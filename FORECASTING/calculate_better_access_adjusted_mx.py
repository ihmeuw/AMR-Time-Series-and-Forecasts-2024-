"""Script to calculate better access scenario based on passed arguments."""

import os
import sys
from functools import partial
from pathlib import Path
from typing import Generator, List

import click
import numpy as np
import xarray as xr
from amr_lib.amr_utilities import args_from_ctx
from amr_lib.fhs_console_utils import run_daly_console, run_yll_console
from db_queries import get_cause_metadata, get_location_metadata
from fhs_lib_cli_tools.lib.fhs_click_decorators import with_config_file
from fhs_lib_data_aggregation.lib.aggregation_methods import summation
from fhs_lib_data_aggregation.lib.aggregator import Aggregator
from fhs_lib_data_aggregation.lib.cause_hierarchy import get_cause_hierarchy
from fhs_lib_data_transformation.lib.dimension_transformation import expand_dimensions
from fhs_lib_database_interface.lib.constants import AgeConstants
from fhs_lib_file_interface.lib.file_interface import FBDPath
from fhs_lib_file_interface.lib.file_system_manager import FileSystemManager
from fhs_lib_file_interface.lib.os_file_system import OSFileSystem
from fhs_lib_file_interface.lib.version_metadata import VersionMetadata
from fhs_lib_file_interface.lib.versioning import Versions
from fhs_lib_file_interface.lib.xarray_wrapper import open_xr, save_xr
from fhs_lib_orchestration_interface.lib.slicer import Slicer
from fhs_lib_orchestration_interface.lib.submit_cluster_job import (
    submit_cluster_job_request,
)
from fhs_lib_year_range_manager.lib import YearRange
from tiny_structured_logger.lib import fhs_logging

FileSystemManager.set_file_system(OSFileSystem())

LOCATION_QUERY = "level==3"  # Query used to subset locations
MX_SUMMARY_FOLDER = "summary"  # Also could be summary_agg

CFR_AGE_TO_DETAILED_AGE_DICT = {
    25: [15, 16, 17, 18],
    26: [19, 20, 30, 31, 32, 235],
    42: [2, 3],
    179: [4, 5],
    239: list(range(6, 15)),
    1: [2, 3, 4, 5],
}


FRAC_SYN_TO_CFR_SYN_DICT = {
    "L2_lower_respiratory_infection_hosp": "lower_respiratory_infection_hospital",
    "L2_encephalitis": "MMO_encephalitis_all",
    "L2_myelitis_meningoencephalitis_other": "MMO_encephalitis_all",
    "L1_bone_joint_infection": "bone_joint_infection_all",
    "L2_meningitis": "meningitis_all",
    "L2_genital_infection": "genital_infection_all",
    "L2_eye_infection": "eye_infection_all",
    "L2_endocarditis": "cardiovascular_infection_all",
    "L1.5_myocarditis_pericarditis_carditis": "cardiovascular_infection_all",
    "L1_peritoneal_and_intra_abdomen_infection": "peritoneal_and_intra_abdomen_infection_all",
    "L2_upper_respiratory_infection": "respiratory_non_lower_all",
    "L2_mix_respiratory_infection": "respiratory_non_lower_all",
    "L2_oral_infection": "oral_infection_all",
    "L2_blood_stream_infection": "blood_stream_infection_all",
    "L2_urinary_tract_infection_hosp": "urinary_tract_infection_hospital",
    "L2_skin_infection": "skin_infection_all",
}

COMMUNITY_DICT = {
    "lri": "lower_respiratory_infection_community",
    "urinary_nephritis": "urinary_tract_infection_community",
    "tb": "tb",
}

HIV_TB_ACAUSES = [
    "hiv",
    "tb",
    "lri",
    "urinary_nephritis",
]

FHS_CAUSE_SET_ID = 6
FHS_LOCATION_SET_ID = 39
MD_SEXES = [1, 2]

DEFAULT_LEX_VERSION = "FILEPATH:0"
DEFAULT_YLD_VERSION = "FILEPATH:0"
DEFAULT_LEX_GBD_ROUND_ID = 7


@click.group()
@with_config_file(str(Path(os.curdir) / "configs"))
@click.option(
    "--rr-version",
    type=FBDPath,
    required=True,
    help="Path to relative risks output by interpolate_better_access_rr.py",
)
@click.option(
    "--population-version",
    type=VersionMetadata.default_parser(
        default_epoch="future", default_stage="population"
    ),
    required=True,
    help="Population version to use for aggregation",
)
@click.option(
    "--syndrome-fraction-version",
    type=VersionMetadata.default_parser(default_epoch="past", default_stage="death"),
    required=True,
    help="Syndrome fraction version output by ETL code",
)
@click.option(
    "--mx-version",
    type=VersionMetadata.default_parser(default_epoch="future", default_stage="death"),
    required=True,
    help="Mortality S8 rates version to apply better access scenario to",
)
@click.option(
    "--output-version",
    type=VersionMetadata.default_parser(default_epoch="future", default_stage="death"),
    required=True,
    help="Version to output better access version to.",
)
@click.option("--draws", type=int, default=500, help="Number of draws in output file")
@click.option(
    "--output-scenario",
    type=int,
    required=True,
    help="Scenario to use for output files",
)
@click.option(
    "--years",
    type=YearRange.parse_year_range,
    default="1990:2021:2050",
    help="Years format",
)
@click.option(
    "--target-year",
    type=int,
    default=2030,
    help="Target year for better access scenario",
)
@click.option("--release-id", type=int, help="release ID for aggregator to use")
@click.option("--gbd-round-id", type=int, help="gbd_round_id for aggregator to use")
@click.option("--queue", type=click.Choice(["all", "long"]), default="all")
@click.option("--acause", type=str, default=None)
@click.option(
    "--yll-daly",
    is_flag=True,
    default=False,
    help="Flag to calculate ylls and dalys from adjusted mx",
)
@click.pass_context
def cli(
    ctx: click.Context,
    queue: str,
    acause: str,
    rr_version: FBDPath,
    population_version: VersionMetadata,
    syndrome_fraction_version: VersionMetadata,
    mx_version: VersionMetadata,
    output_version: VersionMetadata,
    years: YearRange,
    target_year: int,
    release_id: int,
    gbd_round_id: int,
    output_scenario: int,
    yll_daly: bool,
    draws: int,
) -> None:
    """Initialize main CLI function to parse args and pass them to the subcommands.

    Args:
        ctx (click.Context): Context object passed from click
        queue (str): Queue to submit jobs to
        acause (str): Cause to calculate better access scenario for
        rr_version (FBDPath): Path to relative risks output by interpolate_better_access_rr.py
        population_version (VersionMetadata): Population version to use for aggregation
        syndrome_fraction_version (VersionMetadata): Syndrome fraction version output by ETL code
        mx_version (VersionMetadata): Mortality S8 rates version to apply better access scenario to
        output_version (VersionMetadata): Version to output better access version to.
        years (YearRange): Years format
        target_year (int): Target year for better access scenario
        release_id (int): release ID for aggregator to use
        gbd_round_id (int): gbd_round_id for aggregator to use
        output_scenario (int): Scenario to use for output files
        yll_daly (bool): Flag to calculate ylls and dalys from adjusted mx
        draws (int): Number of draws in output file

    Returns:
        None
    """
    if output_version.scenario is None:
        output_version = output_version.with_scenario(output_scenario)
    ctx.obj = {
        "gbd_round_id": gbd_round_id,
        "release_id": release_id,
        "target_year": target_year,
        "years": years,
        "output_version": output_version.default_data_source(gbd_round_id),
        "output_scenario": output_scenario,
        "mx_version": mx_version.default_data_source(gbd_round_id),
        "syndrome_fraction_version": syndrome_fraction_version.default_data_source(
            gbd_round_id
        ),
        "population_version": population_version.default_data_source(gbd_round_id),
        "rr_version": rr_version,
        "queue": queue,
        "draw_slicer": Slicer(years=years, draw=list(range(0, draws))),
        "acause": acause,
        "yll_daly": yll_daly,
        "draws": draws,
    }


@cli.command()
@click.pass_context
def run_yll_daly(ctx: click.Context) -> None:
    """Run YLL and DALY console."""
    logger = fhs_logging.get_logger()
    logger.info("Entering `parallelize_main` function.")

    script = Path(__file__)

    queue = ctx.obj["queue"]
    gbd_round_id = ctx.obj["gbd_round_id"]
    release_id = ctx.obj["release_id"]
    target_year = ctx.obj["target_year"]
    years = ctx.obj["years"]
    output_version = ctx.obj["output_version"]
    mx_version = ctx.obj["mx_version"]
    syndrome_fraction_version = ctx.obj["syndrome_fraction_version"]
    population_version = ctx.obj["population_version"]
    rr_version = ctx.obj["rr_version"]
    output_scenario = ctx.obj["output_scenario"]
    yll_daly = ctx.obj["yll_daly"]
    draws = ctx.obj["draws"]
    yll_jid = run_yll_console(
        versions=Versions(
            str(output_version),
            str(population_version),
            DEFAULT_LEX_VERSION,
            str(output_version.with_stage("yll")),
        ),
        output_scenario=output_scenario,
        gbd_round_id=gbd_round_id,
        draws=draws,
        reference_lex_gbd_round_id=DEFAULT_LEX_GBD_ROUND_ID,
        years=years,
        holds=[1],
        memory=10,
        job_name="better_access_calc_ylls",
    )

    daly_jid = run_daly_console(
        versions=Versions(
            str(output_version.with_stage("yll")),
            str(output_version.with_stage("daly")),
            str(population_version),
            DEFAULT_YLD_VERSION,
        ),
        output_scenario=output_scenario,
        gbd_round_id=gbd_round_id,
        years=years,
        draws=draws,
        holds=[f"afterok:{yll_jid}"],
        extra_flags=["--kill-on-invalid-dep=yes"],
        memory=10,
        job_name="better_access_calc_dalys",
    )


@cli.command()
@click.pass_context
def parallelize_main(ctx: click.Context) -> int:
    """Calculate better access scenario in parallel by draw."""
    logger = fhs_logging.get_logger()
    logger.info("Entering `parallelize_main` function.")

    script = Path(__file__)

    queue = ctx.obj["queue"]
    gbd_round_id = ctx.obj["gbd_round_id"]
    release_id = ctx.obj["release_id"]
    target_year = ctx.obj["target_year"]
    years = ctx.obj["years"]
    output_version = ctx.obj["output_version"]
    mx_version = ctx.obj["mx_version"]
    syndrome_fraction_version = ctx.obj["syndrome_fraction_version"]
    population_version = ctx.obj["population_version"]
    rr_version = ctx.obj["rr_version"]
    output_scenario = ctx.obj["output_scenario"]
    yll_daly = ctx.obj["yll_daly"]
    draws = ctx.obj["draws"]

    compute_memory = 150
    comp_args = args_from_ctx(ctx.obj)
    comp_better_access_command = f"{sys.executable} {script} {comp_args} compute-main"

    comp_better_access_jid = submit_cluster_job_request(
        job_name="compute_better_access",
        cores=1,
        memory=compute_memory,
        runtime="1:00:00",
        cluster_queue=f"{queue}.q",
        number_of_array_tasks=draws,
        command_to_run=comp_better_access_command,
    )

    merge_memory_per_draw = 0.35
    merge_memory = int(np.ceil(merge_memory_per_draw * draws))
    mx_acauses = _get_acauses_from_mx_summary_file(mx_version.data_path())
    acause_tasks = []
    for acause in mx_acauses:
        ctx.obj["acause"] = acause
        acause_args = args_from_ctx(ctx.obj)
        merge_command = f"{sys.executable} {script} " f"{acause_args} " f"merge-main"
        jid = submit_cluster_job_request(
            job_name=f"merge_{acause}",
            job_holds=[comp_better_access_jid],
            cores=1,
            memory=merge_memory,
            runtime="5:00:00",
            cluster_queue=f"{queue}.q",
            command_to_run=merge_command,
        )
        acause_tasks.append(jid)

    if yll_daly:
        yll_jid = run_yll_console(
            versions=Versions(
                str(output_version),
                str(population_version),
                DEFAULT_LEX_VERSION,
                str(output_version.with_stage("yll")),
            ),
            output_scenario=output_scenario,
            gbd_round_id=gbd_round_id,
            draws=draws,
            reference_lex_gbd_round_id=DEFAULT_LEX_GBD_ROUND_ID,
            years=years,
            holds=["afterok:" + ",".join(map(str, acause_tasks))],
            extra_flags=["--kill-on-invalid-dep=yes"],
            memory=10,
            job_name="better_access_calc_ylls",
        )

        daly_jid = run_daly_console(
            versions=Versions(
                str(output_version.with_stage("yll")),
                str(output_version.with_stage("daly")),
                str(population_version),
                DEFAULT_YLD_VERSION,
            ),
            output_scenario=output_scenario,
            gbd_round_id=gbd_round_id,
            years=years,
            draws=draws,
            holds=[f"afterok:{yll_jid}"],
            extra_flags=["--kill-on-invalid-dep=yes"],
            memory=10,
            job_name="better_access_calc_dalys",
        )


@cli.command()
@click.pass_context
def compute_main(ctx: click.Context) -> None:
    """Compute better access scenario.

    Args:
        ctx (click.Context): ctx object passed from cli()
    """
    gbd_round_id = ctx.obj["gbd_round_id"]
    release_id = ctx.obj["release_id"]
    target_year = ctx.obj["target_year"]
    years = ctx.obj["years"]
    output_version = ctx.obj["output_version"]
    mx_version = ctx.obj["mx_version"]
    syndrome_fraction_version = ctx.obj["syndrome_fraction_version"]
    population_version = ctx.obj["population_version"]
    rr_version = ctx.obj["rr_version"]
    output_scenario = ctx.obj["output_scenario"]
    draw = ctx.obj["draw_slicer"].slice()["draw"]
    draws = ctx.obj["draws"]

    cause_metadata = get_cause_metadata(
        release_id=release_id, cause_set_id=FHS_CAUSE_SET_ID
    )
    most_detailed_causes = (
        cause_metadata.query("most_detailed == 1")["acause"].tolist() + HIV_TB_ACAUSES
    )

    (
        relative_reductions,
        cause_syndrome_fractions,
        cause_specific_mx,
        population,
    ) = read_input_data(ctx, draw, most_detailed_causes)

    cause_age_specific_reductions = convert_syndrome_reductions_to_cause_reductions(
        ctx, relative_reductions, cause_syndrome_fractions
    )

    # Clip to 1 if any values > 1
    if (cause_age_specific_reductions > 1).any():
        print("Reductions greater than 1. Clipping to 1")
        cause_age_specific_reductions = cause_age_specific_reductions.clip(max=1)

    print("Calculating better access scenario mx.")
    adjusted_mx = cause_specific_mx.sel(location_id=relative_reductions.location_id) * (
        1 - cause_age_specific_reductions
    )

    # Causes from reference without syndrome fraction data
    missing_causes = [
        x
        for x in most_detailed_causes
        if (x not in adjusted_mx["acause"]) and (x in cause_specific_mx["acause"])
    ]

    missing_causes_da = cause_specific_mx.sel(acause=missing_causes)

    adj_and_reference_mx = xr.concat([adjusted_mx, missing_causes_da], dim="acause")

    most_detailed_causes = adj_and_reference_mx.acause.values

    cause_hierarchy = get_cause_hierarchy(most_detailed_causes, gbd_round_id)

    sum_causes = partial(summation, "acause")

    agged_causes = Aggregator.aggregate_entities(
        "acause", cause_hierarchy, adj_and_reference_mx, func=sum_causes
    )

    agged_causes = expand_dimensions(agged_causes, scenario=[output_scenario])

    # Demographic aggregation
    demographic_aggregates = Aggregator.aggregate_everything(
        pop=population, gbd_round_id=gbd_round_id, data=agged_causes
    )
    agg_rate = demographic_aggregates.data.rate
    agg_num = demographic_aggregates.data.number

    out_path_rate = output_version.data_path() / "draw_files" / f"{draw}.nc"
    save_xr(agg_rate, out_path_rate, metric="rate", space="identity")
    print(f"Data saved to {out_path_rate}")

    out_path_num = (
        output_version.append_version_suffix("_num").data_path()
        / "draw_files"
        / f"{draw}.nc"
    )
    save_xr(agg_num, out_path_num, metric="number", space="identity")
    print(f"Data saved to {out_path_num}")

    print("Done!")


@cli.command()
@click.pass_context
def merge_main(ctx: click.Context) -> None:
    """Merge computed values by draw.

    Args:
        ctx (click.Context): ctx object passed from cli()
    """
    gbd_round_id = ctx.obj["gbd_round_id"]
    release_id = ctx.obj["release_id"]
    target_year = ctx.obj["target_year"]
    years = ctx.obj["years"]
    output_version = ctx.obj["output_version"]
    output_scenario = ctx.obj["output_scenario"]
    mx_version = ctx.obj["mx_version"]
    syndrome_fraction_version = ctx.obj["syndrome_fraction_version"]
    population_version = ctx.obj["population_version"]
    rr_version = ctx.obj["rr_version"]
    acause = ctx.obj["acause"]
    draws = range(ctx.obj["draws"])

    rate_da = _load_access_mortality_draws(acause, draws, output_version.data_path())
    save_xr(
        rate_da,
        output_version.data_dir() / f"{acause}.nc",
        metric="rate",
        space="identity",
    )
    print(f"{acause} data saved to {output_version.data_path}")

    del rate_da

    num_da = _load_access_mortality_draws(
        acause, draws, output_version.append_version_suffix("_num").data_path()
    )
    save_xr(
        num_da,
        output_version.append_version_suffix("_num").data_path() / f"{acause}.nc",
        metric="number",
        space="identity",
    )
    print(
        f"{acause} data saved to {output_version.append_version_suffix('_num').data_path()}"
    )


def read_input_data(
    ctx: click.Context, draw: int, most_detailed_causes: List[int]
) -> tuple:
    """Read all input data in.

    Args:
        ctx (click.Context): ctx object passed from cli()
        draw (int): draw to read in
        most_detailed_causes (list): list of most detailed causes to read in

    Returns:
        tuple: Tuple of 4 DataArrays - relative_reductions, cause_syndrome_fractions, cause_specific_mx, population
    """
    gbd_round_id = ctx.obj["gbd_round_id"]
    release_id = ctx.obj["release_id"]
    target_year = ctx.obj["target_year"]
    years = ctx.obj["years"]
    output_version = ctx.obj["output_version"]
    output_scenario = ctx.obj["output_scenario"]
    mx_version = ctx.obj["mx_version"]
    syndrome_fraction_version = ctx.obj["syndrome_fraction_version"]
    population_version = ctx.obj["population_version"]
    rr_version = ctx.obj["rr_version"]
    acause = ctx.obj["acause"]

    locs = (
        get_location_metadata(
            release_id=release_id, location_set_id=FHS_LOCATION_SET_ID
        )
        .query(LOCATION_QUERY)
        .location_id.values
    )

    relative_reductions = (
        xr.open_mfdataset(
            rr_version.glob("*.nc"),
            concat_dim="syndrome",
            combine="nested",
        )["cfr"]
        .sel(year_id=target_year, location_id=locs, drop=True)
        .load()
    )

    cause_syndrome_fractions = (
        xr.open_mfdataset(
            syndrome_fraction_version.data_path().glob("*.nc"),
            concat_dim="acause",
            combine="nested",
            preprocess=_merge_preprocess(years.past_end)._select_past_end_drop_cause_id,
        )["value"]
        .sel(
            age_group_id=AgeConstants.MOST_DETAILED_AGE_GROUP_IDS_PRE_GBD_2020,
            sex_id=MD_SEXES,
            location_id=locs,
        )
        .load()
    )

    relative_reductions = relative_reductions.sel(
        syndrome=list(
            set(list(FRAC_SYN_TO_CFR_SYN_DICT.values()) + list(COMMUNITY_DICT.values()))
        ),
    )
    cause_syndrome_fractions = cause_syndrome_fractions.sel(
        syndrome=list(FRAC_SYN_TO_CFR_SYN_DICT.keys()),
        location_id=relative_reductions.location_id,
    )

    # clip to 1
    cause_syndrome_fractions = cause_syndrome_fractions.clip(max=1)

    cause_specific_mx = _load_reference_mx_data(
        draw, most_detailed_causes, mx_version, output_scenario, locs
    )

    population = (
        open_xr(population_version.data_path() / "population.nc")
        .sel(
            scenario=(
                population_version.scenario
                if population_version.scenario is not None
                else output_scenario
            ),
            location_id=locs,
            drop=True,
        )
        .sel(draw=[draw])
    )

    return relative_reductions, cause_syndrome_fractions, cause_specific_mx, population


def convert_syndrome_reductions_to_cause_reductions(
    ctx: click.Context,
    relative_reductions: xr.DataArray,
    cause_syndrome_fractions: xr.DataArray,
) -> xr.DataArray:
    """Convert the relative reductions to by cause instead of by syndrome.

    Args:
        ctx (click.Context): ctx object passed from cli()
        relative_reductions (xr.DataArray): DataArray of relative reductions by syndrome.
        cause_syndrome_fractions (xr.DataArray): DataArray of proportion of syndromes attributable to a cause.

    Returns:
        xr.DataArray: Relative Reductions by cause.
    """
    gbd_round_id = ctx.obj["gbd_round_id"]
    release_id = ctx.obj["release_id"]
    target_year = ctx.obj["target_year"]
    years = ctx.obj["years"]
    output_version = ctx.obj["output_version"]
    output_scenario = ctx.obj["output_scenario"]
    mx_version = ctx.obj["mx_version"]
    syndrome_fraction_version = ctx.obj["syndrome_fraction_version"]
    population_version = ctx.obj["population_version"]
    rr_version = ctx.obj["rr_version"]
    acause = ctx.obj["acause"]

    syndrome_fractions = []
    logger = fhs_logging.get_logger()

    for fraction_syndrome, cfr_syndrome in FRAC_SYN_TO_CFR_SYN_DICT.items():
        print(f"{fraction_syndrome} * {cfr_syndrome}")

        syndrome_relative_reduction = relative_reductions.sel(
            syndrome=cfr_syndrome, drop=True
        ).dropna("age_group_id", how="all")

        age_fractions = []
        for age_group in syndrome_relative_reduction.age_group_id.values:
            detailed_ages = CFR_AGE_TO_DETAILED_AGE_DICT[age_group]

            try:
                syndrome_specific_fractions = cause_syndrome_fractions.sel(
                    syndrome=fraction_syndrome
                ).fillna(0)

                detailed_age_intersection = list(
                    set(AgeConstants.MOST_DETAILED_AGE_GROUP_IDS_PRE_GBD_2020)
                    & set(detailed_ages)
                )

                reduced_fraction = (
                    syndrome_relative_reduction.sel(age_group_id=age_group, drop=True)
                    * syndrome_specific_fractions.sel(
                        age_group_id=detailed_age_intersection
                    )
                ).persist()

                age_fractions.append(reduced_fraction)
            except Exception as e:
                logger.warn(
                    f"Age group {age_group} does not exist for {fraction_syndrome}, skipping"
                )
                continue

        syndrome_fraction = xr.concat(age_fractions, dim="age_group_id")
        syndrome_fractions.append(syndrome_fraction)

    for cause, cfr_syndrome in COMMUNITY_DICT.items():
        print(f"{cause} * {cfr_syndrome}")

        syndrome_relative_reduction = relative_reductions.sel(
            syndrome=cfr_syndrome, drop=True
        ).dropna("age_group_id", how="all")
        age_fractions = []
        for age_group in syndrome_relative_reduction.age_group_id.values:
            detailed_ages = CFR_AGE_TO_DETAILED_AGE_DICT[age_group]

            try:
                syndrome_specific_fractions = expand_dimensions(
                    xr.zeros_like(
                        cause_syndrome_fractions.sel(
                            acause=cause_syndrome_fractions.sel(
                                syndrome=list(FRAC_SYN_TO_CFR_SYN_DICT.keys())[0]
                            ).acause.values.tolist()[0],
                            syndrome=list(FRAC_SYN_TO_CFR_SYN_DICT.keys())[0],
                            drop=True,
                        ),
                    ).where(False, 1),
                    acause=[cause],
                    syndrome=[cfr_syndrome],
                ).sel(syndrome=cfr_syndrome)
                detailed_age_intersection = list(
                    set(AgeConstants.MOST_DETAILED_AGE_GROUP_IDS_PRE_GBD_2020)
                    & set(detailed_ages)
                )
                reduced_fraction = (
                    syndrome_relative_reduction.sel(age_group_id=age_group, drop=True)
                    * syndrome_specific_fractions.sel(
                        age_group_id=detailed_age_intersection
                    )
                ).persist()
                age_fractions.append(reduced_fraction)
            except Exception as e:
                logger.warn(
                    f"Age group {age_group} does not exist for {fraction_syndrome}, skipping. Error : \n {e}"
                )
                continue

        syndrome_fraction = xr.concat(age_fractions, dim="age_group_id")
        syndrome_fractions.append(syndrome_fraction)

    reduction_target = (
        xr.concat(syndrome_fractions, dim="syndrome")
        .sum("syndrome")
        .load()
        .assign_coords(year_id=2030)
    )
    no_reduction = xr.zeros_like(reduction_target).assign_coords(year_id=2021)

    reduction_interpolated = xr.concat(
        [no_reduction, reduction_target], dim="year_id"
    ).interp(year_id=range(years.past_end + 1, target_year), method="linear")

    reduction_post_target = expand_dimensions(
        reduction_target, year_id=list(range(target_year + 1, years.forecast_end + 1))
    )

    reduction_timeseries = xr.concat(
        [
            reduction_interpolated,
            reduction_target.expand_dims("year_id"),
            reduction_post_target,
        ],
        dim="year_id",
    ).rename("value")

    return reduction_timeseries


class _merge_preprocess:
    def __init__(self, past_end: int) -> None:
        self.past_end = past_end

    def _select_past_end_drop_cause_id(self, da: xr.DataArray) -> xr.DataArray:
        return da.sel(year_id=self.past_end, cause_id=da.cause_id.item(), drop=True)


def _get_acauses_from_mx_summary_file(mx_path: FBDPath) -> np.array:
    return open_xr(mx_path / "summary.nc").acause.values


def _load_reference_mx_data(
    draw: int,
    most_detailed_causes: list,
    mx_version: VersionMetadata,
    output_scenario_id: int,
    locs: List[int],
) -> xr.DataArray:
    all_mx_acauses = [path.stem for path in mx_version.data_path().glob("*.nc")]
    most_detailed_mx_causes = set(all_mx_acauses) & set(most_detailed_causes)

    acause_das = []
    for acause in most_detailed_mx_causes:
        print(acause)

        da = (
            open_xr(mx_version.data_path() / f"{acause}.nc")
            .sel(
                draw=[draw],
                scenario=(
                    mx_version.scenario
                    if mx_version.scenario is not None
                    else output_scenario_id
                ),
                location_id=locs,
                age_group_id=AgeConstants.MOST_DETAILED_AGE_GROUP_IDS_PRE_GBD_2020,
                sex_id=MD_SEXES,
                drop=True,
            )
            .assign_coords(acause=acause)
        )

        acause_das.append(da)

    print("Concatenating acause data.")
    one_draw_da = xr.concat(acause_das, dim="acause", fill_value=0)

    return one_draw_da


def _load_access_mortality_draws(
    acause: str, draws: Generator, path: FBDPath
) -> xr.DataArray:
    logger = fhs_logging.get_logger()

    draw_das = []
    for draw in draws:
        logger.info(f"draw {draw}")

        da = open_xr(path / "draw_files" / f"{draw}.nc")
        if acause not in da.acause.values:
            raise ValueError(f"{acause} not found in draw files.")
        else:
            da = da.sel(acause=acause, drop=True)
        draw_das.append(da)

    print(f"Concatenating {acause} draws from {path}.")
    one_cause_da = xr.concat(draw_das, dim="draw")

    return one_cause_da


if __name__ == "__main__":
    cli()
