"""Code to run attributable burden for AMR.

NOTE:
    - parallel-compute-attributable-burdens ignores save
    version and saves as YYYYMMDD_attributable_burden_STAGE_PAF_VERSION
    - If no --save-version is passed to compute-attributable-burden
    the version name will automatically be set to YYYYMMDD_amr_attributable_STAGE
    where YYYYMMDD is todays date and STAGE is the passed in stage.

example calls:
bash::
# Run mean-level attributable burden code in parallel
# for Death, DALY, YLL, and YLD with
# DALY depending on YLL and YLD attributable
# Burden output
python calculate_attributable_burden.py \
    --gbd-round-id 7 \
    --release-id 9 \
    --versions FILEPATH \
    -v FILEPATH \
    -v FILEPATH \
    -v FILEPATH \
    --stage death \
    --stage yll \
    --stage yld \
    --stage daly \
    --epoch future \
    parallel-compute-attributable-burdens

# Run draw-level attributable burden code in parallel
# for Death, DALY, YLL, and YLD with
# DALY depending on YLL and YLD attributable
# Burden output
python calculate_attributable_burden.py \
    --gbd-round-id 7 \
    --release-id 9 \
    --versions FILEPATH \
    -v FILEPATH  \
    -v FILEPATH \
    -v FILEPATH \
    --stage death \
    --stage yll \
    --stage yld \
    --stage daly \
    --epoch future \
    --draws \
    --num-draws 500 \
    parallel-compute-attributable-burdens


# Get Attributable Burden for all AMR-Cause Risk Pairs

Mean Level DEATH:
python attributable_burden_amr.py \
    --gbd-round-id 7 \
    --release-id 9 \
    --versions FILEPATH \
    -v FILEPATH \
    -v FILEPATH \
    --cur-stage death \
    --epoch future \
    --save-version VERSION \
    compute-attributable-burden


Draw Level DEATH:
python attributable_burden_amr.py \
    --gbd-round-id 7 \
    --release-id 9 \
    --versions FILEPATH \
    -v FILEPATH \
    -v FILEPATH \
    --cur-stage death \
    --epoch future \
    --save-version VERSION \
    --draws \
    --num-draws 500 \
    compute-attributable-burden



# Get Attributable Burden for one AMR-Cause Risk Pair
# NOTE must be in the PAF directory as acause.nc

python calculate_attributable_burden.py \
    --gbd-round-id 7 \
    --release-id 9 \
    --versions FILEPATH \
    -v FILEPATH \
    -v FILEPATH \
    --cur-stage yll \
    --epoch future \
    --save-version VERSION \
    --draws \
    --num-draws 500 \
    --acause ACAUSE \
    calculate-one-cause-main

"""

import sys
from functools import partial
from pathlib import Path
from typing import List, Optional, Union

import click
import db_queries
import xarray as xr
from amr_lib.amr_utilities import (
    args_from_ctx,
    ctx_tuple_to_list,
    date_now,
    filter_dataarrays,
    filter_scenario,
    get_mean_level_pafs,
    get_most_detailed,
    get_one_paf,
)
from amr_lib.fhs_console_utils import run_aggregator_console, run_summary_maker_console
from fhs_lib_data_aggregation.lib import aggregator
from fhs_lib_data_aggregation.lib.aggregation_methods import summation
from fhs_lib_data_aggregation.lib.aggregator import Aggregator
from fhs_lib_data_aggregation.lib.cause_hierarchy import get_cause_hierarchy
from fhs_lib_data_transformation.lib.resample import resample
from fhs_lib_file_interface.lib.file_interface import FBDPath
from fhs_lib_file_interface.lib.file_system_manager import FileSystemManager
from fhs_lib_file_interface.lib.os_file_system import OSFileSystem
from fhs_lib_file_interface.lib.versioning import Versions
from fhs_lib_file_interface.lib.xarray_wrapper import open_xr, save_xr
from fhs_lib_orchestration_interface.lib.submit_cluster_job import (
    submit_cluster_job_request,
)
from fhs_lib_year_range_manager.lib.year_range import YearRange
from tiny_structured_logger.lib import fhs_logging

FileSystemManager.set_file_system(OSFileSystem())

logger = fhs_logging.get_logger()
REFERENCE_SCENARIO = 0
FHS_CAUSE_SET = 6
FHS_LOCATION_SET = 39
EPOCHS = ["past", "future"]
BURDENS = ["death", "yll", "yld", "daly"]
DALY_VERSION_DEPENDENCIES = ["population", "yll", "yld"]
TOTALS_REQUIRED = ["yll", "death"]
VERSION_DEPENDENCIES = ["population", "paf"]
REPORTING_CAUSES = "level <= 2"

JOB_MEM = 50
JOB_TIME = "3:00:00"
JOB_QUEUE = "all.q"
JOB_CORES = 10

YLL = "yll"
YLD = "yld"
YEARS = YearRange.parse_year_range("1990:2022:2050")

R7_MOST_DETAILED_NOT_IN_R6 = {
    "endo_other",
    "digest_ibd_colitis",
    "endo_thyroid",
    "inj_electrocution",
    "digest_ibd_crohns",
    "neo_ben_other",
    "mental_eating_anorexia",
}

ATTRS_RATE = {"metric": "rate", "space": "identity"}
ATTRS_NUM = {"metric": "rate", "space": "identity"}

ATTRIBUTABLE = "attributable"
ASSOCIATED = "associated"

# we want to remove lri_corona from our burdens since we're not modeling covid
RM_CHILDREN = {"_ri": "lri_corona"}

AGGREGATOR_DEFAULTS = {
    "log_level": "DEBUG",
    "provenance": False,
    "save_number": True,
    "save_percent": False,
    "cause_aggregation": True,
    "demographic_aggregation": True,
    "rake_subnats": False,
    "national_only": True,
    "custom_cause_hierarchies": [],
    "custom_location_set_ids": [],
}

SUMMARYMAKER_DEFAULTS = {
    "log_level": "DEBUG",
    "national_only": False,
    "combine_past_and_future": False,
    "compute_demographic_aggregates": False,
    "entity_type": "acause",
    "input_entities_aggregated": True,
    "no_missing_entity_check": True,
}


def ensure_rate_space(
    burden_da: xr.DataArray, arg_object: Union[click.Context, dict]
) -> None:
    """Ensure that the passed dataArray is in rate space.

    Args:
        burden_da (xr.DataArray): DataArray to ensure is in rate space
        arg_object (Union[click.Context, dict]):  object of args being passed from cli()

    """
    if burden_da.metric != "rate":
        raise ValueError(
            (
                "Burden Version Passed Must be in rate space!!!"
                f"Version {arg_object['burden_version']} for "
                f"{arg_object['stage']} in {burden_da.metric} space."
            )
        )


def ensure_correct_versions(ctx: dict) -> None:
    """Make sure that the correct versions are given for called stage.

    Args:
        ctx (dict): ctx object created in cli()

    Returns:
        None

    Raises:
        ValueError - when missing needed version
    """
    versions = ctx["versions"]
    stage = ctx["cur_stage"]
    epoch = ctx["epoch"]
    if stage is None:
        for x in ctx["stage"]:
            arg = ctx.copy()
            arg["cur_stage"] = x
            ensure_correct_versions(arg)
            return
    # two cases - dalys and not dalys
    if stage == "daly":
        if (
            sum(
                [
                    x
                    not in list(set(list(versions[epoch].keys()) + list(ctx["stage"])))
                    for x in DALY_VERSION_DEPENDENCIES
                ]
            )
            > 0
        ):
            raise ValueError(
                "Must include version for population, yll, and yld for DALYs stage"
            )
    elif stage == "yld":
        if (
            sum(
                [
                    x
                    not in list(set(list(versions[epoch].keys()) + list(ctx["stage"])))
                    for x in ["population", "yll"]
                ]
            )
            > 0
        ):
            raise ValueError(
                f"Must include version for population, paf, and yll for {stage}'s stage"
            )
    elif stage == "all":
        if (
            sum(
                [
                    x not in versions[epoch].keys()
                    for x in TOTALS_REQUIRED + VERSION_DEPENDENCIES
                ]
            )
            > 0
        ):
            raise ValueError(
                "Must include version for population, yll, and yld for DALYs stage"
            )
    else:
        if (
            sum(
                [
                    x not in versions[epoch].keys()
                    for x in VERSION_DEPENDENCIES + [stage]
                ]
            )
            > 0
        ):
            raise ValueError(
                f"Must include version for population, paf, and {stage} for {stage}'s stage"
            )


def make_save_version(ctx: click.Context) -> str:
    """Initialize save_version if nothing passed.

    Args:
        ctx (click.context): ctx_object passed by cli()

    returns str

    """
    if ctx.obj["associated"]:
        return (
            f"{date_now()}_amr_associated_{ctx.obj['stage']}{ctx.obj['save_suffix']}"
            if ctx.obj["save_version"] is None
            else f'{ctx.obj["save_version"]}{ctx.obj["save_suffix"]}'
        )
    else:
        return (
            f"{date_now()}_amr_attributable_{ctx.obj['stage']}{ctx.obj['save_suffix']}"
            if ctx.obj["save_version"] is None
            else f'{ctx.obj["save_version"]}{ctx.obj["save_suffix"]}'
        )


def agg_causes(
    da: xr.DataArray,
    release_id: int,
    gbd_round_id: int,
) -> xr.Dataset:
    """Aggregate Causes up the Cause Hierarchy and return reporting causes.

    Args:
        da (xr.DataArray): DataArray with most detailed causes to aggregate
        release_id (int): Release-ID to pull cause-metadata with
        gbd_round_id (int): GBD Round id to get cause-hierarchy from
    """
    all_causes = db_queries.get_cause_metadata(
        release_id=release_id, cause_set_id=FHS_CAUSE_SET
    )

    most_detailed_causes = da.acause.values

    cause_hierarchy = get_cause_hierarchy(most_detailed_causes, gbd_round_id)

    logger.info("Aggregating causes.")
    sum_causes = partial(summation, "acause")
    agged_causes = Aggregator.aggregate_entities(
        "acause", cause_hierarchy, da, func=sum_causes
    )

    reporting_causes = all_causes.query(REPORTING_CAUSES).acause.to_list()
    reporting_causes = list(
        set(reporting_causes).intersection(agged_causes.acause.values)
    )

    return agged_causes.sel(acause=reporting_causes)


def make_summaries(arg_dict: dict, holds: List[int] = [1]) -> None:
    """Make summary maker calls based on passed args.

    Args:
        arg_dict (dict): Dictionary of arguments from cli.ob
        holds (List[int], optional): Jobs to wait to run summary maker for. Defaults to [1].
    """
    gbd_round_id = arg_dict["gbd_round_id"]
    release_id = arg_dict["release_id"]
    scenario = arg_dict["scenario"]
    versions = arg_dict["versions"]
    stage = arg_dict["cur_stage"]
    epoch = arg_dict["epoch"]
    save_version = arg_dict["save_version"]
    draws = arg_dict["draws"]
    associated = arg_dict["associated"]

    args = dict(
        stage=stage,
        gbd_round_id=gbd_round_id,
        years=str(YEARS),
        epoch=epoch,
        holds=holds,
        extra_flags=[] if holds == [1] else ["--kill-on-invalid-dep=yes"],
        queue=f"{JOB_QUEUE.split('.')[0]}",
    )

    run_summary_maker_console(
        job_name=f"summarize_{stage}_draws_nums",
        versions=Versions(str(f"7/{epoch}/{stage}/{save_version}_agg_num")),
        **args,
        **SUMMARYMAKER_DEFAULTS,
    )

    run_summary_maker_console(
        job_name=f"summarize_{stage}_draws_rates",
        versions=Versions(str(f"7/{epoch}/{stage}/{save_version}_agg")),
        **args,
        **SUMMARYMAKER_DEFAULTS,
    )


def compute_attributable_burden_code(
    arg_dict: dict,
    holds: List[int] = [1],
    paf_everywhere: bool = False,
) -> List[int]:
    """Submit jobs to the cluster for every cause for stage specified by arguments.

    Args:
        arg_dict(dict): Dictionary of arguments from cli.obj

    Returns:
        List[int] - List of job ids submitted to the cluster

    """
    gbd_round_id = arg_dict["gbd_round_id"]
    release_id = arg_dict["release_id"]
    scenario = arg_dict["scenario"]
    versions = arg_dict["versions"]
    stage = arg_dict["cur_stage"]
    epoch = arg_dict["epoch"]
    save_version = arg_dict["save_version"]
    draws = arg_dict["draws"]
    associated = arg_dict["associated"]
    round_shift_id = arg_dict["round_shift_id"]

    # add scenario to all versions so
    # aggregator doesn't fail in single scenario mode
    for version in versions:
        if (version.scenario is None) and (epoch != "past"):
            version = version.with_scenario(scenario)
            versions = versions.replace(version.epoch, version.stage, version)

    if not draws:
        raise ValueError("Current Call only works for Draw Level Attributable burden")

    base_args = arg_dict.copy()

    jids = [1]
    command = Path(__file__).absolute()
    if (stage in ["death", "yll"]) or paf_everywhere:
        paf = versions.get_version_metadata(epoch, "paf")
        if not paf.data_dir().exists():
            raise FileNotFoundError(f"PAF directory {paf.data_dir()} does not exist")
        for acause_path in paf.data_dir().glob("*.nc"):
            acause = acause_path.stem
            base_args["save_suffix"] = ""
            base_args["save_version"] = save_version
            base_args["acause"] = acause

            cause_args = args_from_ctx(base_args)
            cause_command = (
                f"{sys.executable} {command} "
                f"{cause_args} "
                " calculate-one-cause-main "
            )

            new_jid = submit_cluster_job_request(
                job_name=f"calculate_{acause}_{'associated' if associated else 'attributable'}_{stage}_draws",
                cores=JOB_CORES,
                memory=int(JOB_MEM),
                runtime=JOB_TIME,
                cluster_queue=JOB_QUEUE,
                command_to_run=cause_command,
                job_holds=holds,
            )

            jids.append(new_jid)

        agg_args = AGGREGATOR_DEFAULTS.copy()
        if epoch == "future":
            agg_args["versions"] = Versions("FILEPATH")
            agg_args["output_version"] = "FILEPATH"
            agg_args["output_scenario"] = scenario
        else:
            agg_args["versions"] = Versions("FILEPATH")
            agg_args["output_version"] = arg_dict["save_version"] + "_agg"
            agg_args["output_scenario"] = scenario
        agg_args["epoch"] = epoch
        agg_args["gbd_round_id"] = (
            round_shift_id if round_shift_id is not None else gbd_round_id
        )
        agg_args["stage"] = stage
        agg_args["years"] = str(YEARS)
        agg_args["draws"] = arg_dict["num_draws"]
        agg_args["queue"] = f"{JOB_QUEUE.split('.')[0]}"

        aggregation_jid = run_aggregator_console(
            **agg_args, job_name=f"aggregate_{stage}_draws", memory=JOB_MEM, holds=jids
        )
        make_summaries(arg_dict, [aggregation_jid])

    else:
        yll = versions.get_version_metadata(epoch, YLL)
        for acause_path in yll.data_dir().glob("*.nc"):
            acause = acause_path.stem
            base_args["save_suffix"] = ""
            base_args["save_version"] = save_version
            base_args["acause"] = acause

            cause_args = args_from_ctx(base_args)
            cause_command = (
                f"{sys.executable} {command} "
                f"{cause_args} "
                " compute-one-cause-main "
            )

            new_jid = submit_cluster_job_request(
                job_name=f"calculate_{acause}_{'associated' if associated else 'attributable'}_{stage}_draws",
                cores=JOB_CORES,
                memory=int(JOB_MEM),
                runtime=JOB_TIME,
                cluster_queue=JOB_QUEUE,
                command_to_run=cause_command,
                job_holds=holds,
            )
            jids.append(new_jid)
        make_summaries(arg_dict, jids)

    return jids


def compute_mean_level_attributable_burden_code(
    arg_dict: dict,
) -> None:
    """Calculate attributable burden at the mean level from arguments.

    Args:
        arg_dict(dict): Dictionary of arguments from cli.obj
    """
    gbd_round_id = arg_dict["gbd_round_id"]
    release_id = arg_dict["release_id"]
    versions = arg_dict["versions"]
    stage = arg_dict["cur_stage"]
    epoch = arg_dict["epoch"]
    save_version = arg_dict["save_version"]
    round_shift_id = arg_dict["round_shift_id"]
    yld_multiplier = arg_dict["yld_multiplier"]

    save_version = FBDPath(f"FILEPATH/{save_version}")

    location_metadata = db_queries.get_location_metadata(
        release_id=release_id, location_set_id=FHS_LOCATION_SET
    )

    nats = location_metadata.query("level == 3").location_id.to_list()
    if stage in ["death", "yll"]:
        # death and YLL take the total deaths/ylls and multiply them
        # by pafs to get attributable burden
        all_pafs = (
            filter_dataarrays(
                arg_dict,
                stage="paf",
                da=get_mean_level_pafs(versions.data_dir(gbd_round_id, epoch, "paf")),
                versions=versions,
            )
            .sel(
                year_id=(
                    YEARS.forecast_years if epoch == "future" else YEARS.past_years
                )
            )
            .sel(location_id=nats)
        )

        burdens = filter_dataarrays(
            arg_dict,
            stage=stage,
            da=open_xr("FILEPATH").chunk("auto"),
            versions=versions,
        )

        ensure_rate_space(burdens, arg_dict)
        logger.info(f"Loading {stage} acauses we have PAFs for ")

        for coord in all_pafs.coords:
            logger.info("Taking out covid since we're not forecasting for AMR")
            coord_value = all_pafs[coord].values.tolist()
            if coord == "acause":
                # taking COVID out since we don't forecast
                for parent, child in RM_CHILDREN.items():
                    burdens = burdens.where(
                        burdens.coords["acause"] != parent,
                        burdens.sel(acause=parent)
                        - burdens.sel(acause=child, drop=True),
                    )
                    coord_value = list(set(coord_value) - R7_MOST_DETAILED_NOT_IN_R6)
            burdens = burdens.sel(**{coord: coord_value})
        logger.info("Calculating")
        attributable_burdens = all_pafs * burdens
    elif stage == "yld":
        # Based on previous results we can estimate that
        # YLDs are 0.5% of YLLs or YLLS * .005
        attributable_burdens = (
            filter_dataarrays(
                arg_dict,
                stage="yll",
                da=open_xr("FILEPATH").chunk("auto"),
                versions=versions,
            )
            * yld_multiplier
        )

    else:
        # Sum YLLs and YLDs to get DALYs
        ylls = filter_dataarrays(
            arg_dict,
            stage="yll",
            da=open_xr("FILEPATH").chunk("auto"),
            versions=versions,
        )
        ylds = filter_dataarrays(
            arg_dict,
            stage="yld",
            da=open_xr("FILEPATH").chunk("auto"),
            versions=versions,
        )

        attributable_burdens = ylls + ylds

    agged_attributable = aggregator.Aggregator.aggregate_everything(
        pop=filter_dataarrays(
            arg_dict,
            stage="population",
            da=open_xr("FILEPATH").chunk("auto"),
            age_filter=False,
            versions=versions,
        ),
        gbd_round_id=round_shift_id if round_shift_id is not None else gbd_round_id,
        data=attributable_burdens,
    )
    agged_rates = agged_attributable.data.rate
    agged_counts = agged_attributable.data.number

    # get aggregated causes for counts and rates
    agged_rates = agg_causes(
        agged_rates,
        (release_id if round_shift_id is not None else 6),
        (round_shift_id if round_shift_id is not None else gbd_round_id),
    )
    agged_counts = agg_causes(
        agged_counts,
        (release_id if round_shift_id is not None else 6),
        (round_shift_id if round_shift_id is not None else gbd_round_id),
    )

    save_xr(agged_rates, (save_version / "summary/summary.nc"), **ATTRS_RATE)
    logger.info(
        f'Aggregated rates saved to  {str(save_version / "summary/summary.nc")}'
    )

    save_xr(
        agged_counts,
        (FBDPath(str(save_version.LFN()) + "_num") / "summary/summary.nc"),
        **ATTRS_NUM,
    )
    logger.info(
        f'Aggregated counts saved to  {str(FBDPath(str(save_version.LFN()) + "_num") / "summary/summary.nc")}'
    )


def compute_one_cause(
    arg_dict: dict,
) -> None:
    """Calculate draw level attributable burden for one cause based on arguments.

    Args:
        arg_dict(dict): Dictionary of arguments from cli.obj

    """
    gbd_round_id = arg_dict["gbd_round_id"]
    release_id = arg_dict["release_id"]
    versions = arg_dict["versions"]
    stage = arg_dict["cur_stage"]
    epoch = arg_dict["epoch"]
    save_version = arg_dict["save_version"]
    acause = arg_dict["acause"]
    num_draws = arg_dict["num_draws"]
    yld_multiplier = arg_dict["yld_multiplier"]

    save_version = FBDPath("FILEPATH")

    if stage in ["death", "yll"]:
        # death and YLL take the total deaths/ylls and multiply them
        # by pafs to get attributable burden
        cause_paf = get_one_paf(
            versions.data_dir(gbd_round_id, epoch, "paf"), acause, draw_level=True
        ).chunk("auto")
        cause_paf = get_most_detailed(
            cause_paf, release_id, nationals_most_detailed=True
        ).sel(year_id=(YEARS.forecast_years if epoch == "future" else YEARS.past_years))
        cause_paf = filter_scenario(arg_dict, "paf", cause_paf, versions)

        burdens = open_xr(
            versions.data_dir(gbd_round_id, epoch, stage) / f"{acause}.nc"
        ).chunk("auto")

        ensure_rate_space(burdens, arg_dict)
        logger.info("Removing un-modeled children")
        if acause in RM_CHILDREN.keys():
            rm_child = open_xr(
                versions.data_dir(gbd_round_id, epoch, stage)
                / f"{RM_CHILDREN[acause]}.nc"
            )
            if "acause" in rm_child.coords:
                rm_child = rm_child.sel(acause=RM_CHILDREN[acause], drop=True)

            burdens = burdens - rm_child

        logger.info("filtering burdens")
        burdens = get_most_detailed(
            burdens, release_id, nationals_most_detailed=True
        ).sel(year_id=(YEARS.forecast_years if epoch == "future" else YEARS.past_years))
        burdens = filter_scenario(arg_dict, stage, burdens, versions)

        burdens = resample(burdens, num_draws)
        cause_paf = resample(cause_paf, num_draws)

        logger.info(f"Loading {stage} acauses we have PAFs for ")

        for coord in cause_paf.coords:
            if coord != "acause":
                coord_value = cause_paf[coord].values.tolist()
                coord_value = list(
                    set(coord_value).intersection(set(burdens[coord].values.tolist()))
                )
                burdens = burdens.sel(**{coord: coord_value})
        logger.info("calculating burden")
        attributable_burdens = cause_paf * burdens
    elif stage == "yld":
        # Based on previous results we can estimate that
        # YLDs are 0.5% of YLLs or YLLS * .005

        yll_da = get_most_detailed(
            open_xr(versions.data_dir(gbd_round_id, epoch, YLL) / f"{acause}.nc"),
            release_id,
            nationals_most_detailed=True,
        )
        yll_da = filter_scenario(arg_dict, YLL, yll_da, versions)
        yll_da = resample(yll_da, num_draws)
        attributable_burdens = yll_da * yld_multiplier

    else:
        # Sum YLLs and YLDs to get DALYs

        yll_da = get_most_detailed(
            open_xr(versions.data_dir(gbd_round_id, epoch, YLL) / f"{acause}.nc"),
            release_id,
            nationals_most_detailed=True,
        )
        yll_da = filter_scenario(arg_dict, YLL, yll_da, versions)
        yll_da = resample(yll_da, num_draws)

        yld_da = get_most_detailed(
            open_xr(versions.data_dir(gbd_round_id, epoch, YLD) / f"{acause}.nc"),
            release_id,
            nationals_most_detailed=True,
        )
        yld_da = filter_scenario(arg_dict, YLD, yld_da, versions)
        yld_da = resample(yld_da, num_draws)

        attributable_burdens = yll_da + yld_da
    logger.info("saving")
    save_xr(attributable_burdens, (save_version / f"{acause}.nc"), **ATTRS_RATE)
    logger.info(f'Un-Aggregated Burdens saved to  {str(save_version / f"{acause}.nc")}')


with_release_id = click.option(
    "--release-id",
    type=int,
    required=False,
    default=9,
    help="Release ID for pulling location sets",
)

with_gbd_round_id = click.option(
    "--gbd-round-id", type=int, required=True, help="GBD Round ID"
)

with_versions = click.option(
    "--versions",
    "-v",
    type=str,
    multiple=True,
    callback=Versions.parse_cli_versions_callback,
    help="Versions to calculate Attributable Burden for.\n"
    "If stage is YLL or Death three versions are needed -  "
    "one for the stages rate values, one for population,"
    "and one for AMR attributable PAFS."
    "If stage is YLD two versions are needed - "
    "Population and Attributable YLL (YLD = 0.005 * YLL)"
    "If stage is DALY three versions are needed - "
    "One for population, One for AMR Attributable YLLs,"
    "and One for AMR Attributable YLDs. \n"
    "These should be in the format 7/future/stage/version "
    ", future/stage/version or stage/version. Note if "
    "stage/version format is used it will be assumed"
    "they are future values.\n"
    "Single scenario specification is enabled"
    "The code assumes the scenario specified as "
    "--scenario, but if single scenario mode selectors"
    "are used for versions that will be used",
)

with_multiplier = click.option(
    "--yld-multiplier",
    required=False,
    default=0.033,
    type=float,
    help="Multiplier to use to compute YLDs from YLLs. Defaults to 0.033",
)

with_stage = click.option(
    "--stage",
    type=str,
    multiple=True,
    required=True,
    help="Burden(s) we are calculating attributable burden for",
)

with_cur_stage = click.option(
    "--cur-stage",
    type=click.Choice(BURDENS),
    required=False,
    help="Burden we are calculating attributable burden for",
)

with_epoch = click.option(
    "--epoch",
    type=click.Choice(EPOCHS),
    required=False,
    default="future",
    help="Epoch we are calculating attributable burden for.\
            Can be past or future.",
)

with_scenario = click.option(
    "--scenario",
    "-s",
    multiple=False,
    default=REFERENCE_SCENARIO,
    help="Scenarios to calculate attributable burden for",
)

with_save_version = click.option(
    "--save-version",
    type=str,
    required=False,
    default=None,
    help="Save version for attributable burden",
)

with_round_shift_id = click.option(
    "--round-shift-id",
    type=int,
    default=None,
    help="GBD Round ID for use round shifting",
)

with_save_suffix = click.option(
    "--save-suffix",
    type=str,
    default="",
    help="Suffix for all save versions for parallel-compute-all-burdens",
)

with_associated = click.option(
    "--associated", is_flag=True, help="Flag to run as attributable burden"
)

with_draws = click.option("--draws", is_flag=True, help="Flag to run with draws")

with_acause = click.option(
    "--acause", type=str, required=False, help="Acause for calculate-one-cause-main"
)

with_num_draws = click.option(
    "--num-draws",
    type=int,
    required=False,
    help="Number of draws to resample to in case of draws. Defaults to 100",
)


@click.group()
@with_release_id
@with_gbd_round_id
@with_versions
@with_stage
@with_cur_stage
@with_epoch
@with_scenario
@with_save_version
@with_save_suffix
@with_associated
@with_acause
@with_draws
@with_num_draws
@with_round_shift_id
@with_multiplier
@click.pass_context
def cli(
    ctx: click.Context,
    versions: Versions,
    stage: str,
    cur_stage: List[str],
    epoch: str,
    gbd_round_id: int,
    release_id: int,
    scenario: List[int],
    associated: bool,
    draws: bool,
    save_version: Optional[FBDPath],
    save_suffix: Optional[str],
    acause: Optional[str],
    num_draws: Optional[int],
    round_shift_id: Optional[int],
    yld_multiplier: Optional[int],
) -> None:
    """Initialize cli object for use in command line calls.

    Args:
        ctx (click.Context): Object passed from CLI interface
        versions (Versions): Versions passed for usage in call
        stage (str): Stage to compute burden for
        epoch (str): Epoch to compute burden for
        gbd_round_id (int): GBD Round ID to use to pull files
        release_id (int): Release ID for db_queries calls
        scenario (List[int]): Scenario to pull as default and set output files as
        associated (bool): Flag to run associated burden instead of attributable
        draws (bool): Flag to run on draws instead of mean level
        save_version (Optional[FBDPath]): Version to save output to
        save_suffix (Optional[str]): Suffix to add to all save versions
        acause (Optional[str]): Acause for use in merge code
        num_draws (Optional[int]): Number of draws to resample files to
        round_shift_id (Optional[int]): GBD Round ID for use in round shifting
        yld_multiplier (Optional[int]): Multiplier to use to compute YLDs from YLLs

    """
    if draws and (num_draws is None):
        num_draws = 100
    ctx.obj = {
        "gbd_round_id": gbd_round_id,
        "release_id": release_id,
        "versions": versions.update_default_data_source(gbd_round_id),
        "stage": stage,
        "cur_stage": (
            cur_stage if ((len(stage) == 1) and (cur_stage is None)) else stage[0]
        ),
        "epoch": epoch,
        "scenario": scenario,
        "save_version": save_version,
        "save_suffix": save_suffix,
        "associated": associated,
        "acause": acause,
        "draws": draws,
        "num_draws": num_draws,
        "round_shift_id": (
            round_shift_id if round_shift_id is not None else gbd_round_id
        ),
        "yld_multiplier": yld_multiplier,
    }

    ctx.obj["save_version"] = make_save_version(ctx)
    ensure_correct_versions(ctx.obj)


@cli.command()
@click.pass_context
def parallel_compute_attributable_burdens(ctx: click.Context) -> None:
    """Compute attributable burdens for all stages in parallel.

    Args:
        ctx (click.Context): ctx object passed from cli()

    """
    ctx = ctx_tuple_to_list(ctx)
    base_args = ctx.obj
    versions = ctx.obj["versions"]
    command = Path(__file__).absolute()
    draws = ctx.obj["draws"]

    # general checks for associated vs attributable
    # to ensure that attributable PAFs aren't used
    # for associated burden and vice versa
    if base_args["associated"]:
        desc_string = ASSOCIATED
        for version in versions:
            if ATTRIBUTABLE in str(version):
                raise ValueError(
                    f"'attributable' should not be in the version name for {str(version)}"
                )
    else:
        desc_string = ATTRIBUTABLE
        for version in versions:
            if ASSOCIATED in str(version):
                raise ValueError(
                    f"'associated' should not be in the version name for {str(version)}"
                )
    stages_to_run = base_args["stage"]

    death_save = (
        f"{date_now()}_{desc_string}_burden_death_"
        f"{versions.get_version_metadata(ctx.obj['epoch'], 'paf').version}{ctx.obj['save_suffix']}"
    )
    yll_save = (
        f"{date_now()}_{desc_string}_burden_yll_",
        f"{versions.get_version_metadata(ctx.obj['epoch'], 'paf').version}",
        f"{ctx.obj['save_suffix']}",
    )
    yld_save = (
        f"{date_now()}_{desc_string}_burden_yld_",
        f"{versions.get_version_metadata(ctx.obj['epoch'], 'paf').version}",
        f"{ctx.obj['save_suffix']}",
    )
    daly_save = (
        f"{date_now()}_{desc_string}_burden_daly_",
        f"{versions.get_version_metadata(ctx.obj['epoch'], 'paf').version}",
        f"{ctx.obj['save_suffix']}",
    )

    base_args["save_suffix"] = ""
    death_args = base_args
    if "death" in stages_to_run:
        death_args["save_version"] = death_save
        death_args["cur_stage"] = "death"

        if draws:
            # Submit job to calculate draw level
            # attributable burden for deaths - parallelized by cause
            _ = compute_attributable_burden_code(death_args, paf_everywhere=True)
        else:
            # Submit job to calculate mean level
            # attributable burden for deaths - 1 job
            death_args = args_from_ctx(death_args)
            death_command = (
                f"{sys.executable} {command} "
                f"{death_args} "
                " compute-attributable-burden "
            )
            _ = submit_cluster_job_request(
                job_name=f"calculate_{desc_string}_death",
                cores=JOB_CORES,
                memory=int(JOB_MEM * 2),
                runtime=JOB_TIME,
                cluster_queue=JOB_QUEUE,
                command_to_run=death_command,
            )

    if "yll" in stages_to_run:
        yll_args = base_args
        yll_args["save_version"] = yll_save
        yll_args["cur_stage"] = "yll"

        if draws:
            # Submit job to calculate draw level
            # attributable burden for YLLs - parallelized by cause
            yll_jobs = compute_attributable_burden_code(death_args, paf_everywhere=True)
        else:
            # Submit job to calculate mean level
            # attributable burden for YLLs - 1 job
            yll_args = args_from_ctx(yll_args)
            yll_command = (
                f"{sys.executable} {command} "
                f"{yll_args} "
                " compute-attributable-burden "
            )

            yll_jid = submit_cluster_job_request(
                job_name=f"calculate_{desc_string}_yll",
                cores=JOB_CORES,
                memory=int(JOB_MEM * 2),
                runtime=JOB_TIME,
                cluster_queue=JOB_QUEUE,
                command_to_run=yll_command,
            )
    else:
        # initialize yll_jobs to 1 if YLLs are not being calculated for job holds
        yll_jobs = [1]

    if "yld" in stages_to_run:
        if "yll" in stages_to_run:
            # if yll run above, make job dependent and pass save version to YLD
            yld_args = base_args
            yld_args["save_version"] = yld_save
            yld_args["cur_stage"] = "yld"
            yld_args["versions"] = Versions(
                str(versions.get_version_metadata(base_args["epoch"], "population")),
                str(versions.get_version_metadata(base_args["epoch"], "paf")),
                f"FILEPATH/{yll_save}",
            )
            if draws:
                yld_jobs = compute_attributable_burden_code(
                    yld_args, yll_jobs, paf_everywhere=True
                )
            else:
                yld_args = args_from_ctx(yld_args)
                yld_command = (
                    f"{sys.executable} {command} "
                    f"{yld_args} "
                    " compute-attributable-burden "
                )

                yld_jid = submit_cluster_job_request(
                    job_name=f"calculate_{desc_string}_yld",
                    cores=JOB_CORES,
                    memory=int(JOB_MEM * 2),
                    runtime=JOB_TIME,
                    cluster_queue=JOB_QUEUE,
                    command_to_run=yld_command,
                    job_holds=[yll_jid],
                )
        elif "yll" in versions[base_args["epoch"]].keys():
            # Otherwise no need to wait and use passed yll version
            yld_args = base_args
            yld_args["save_version"] = yld_save
            yld_args["cur_stage"] = "yld"
            yld_args["versions"] = Versions(
                str(versions.get_version_metadata(base_args["epoch"], "population")),
                str(versions.get_version_metadata(base_args["epoch"], "paf")),
                str(versions.get_version_metadata(base_args["epoch"], "yll")),
            )
            if draws:
                yld_jobs = compute_attributable_burden_code(
                    death_args, [1], paf_everywhere=True
                )
            else:
                yld_args = args_from_ctx(yld_args)
                yld_command = (
                    f"{sys.executable} {command} "
                    f"{yld_args} "
                    " compute-attributable-burden "
                )

                yld_jid = submit_cluster_job_request(
                    job_name=f"calculate_{desc_string}_yld",
                    cores=JOB_CORES,
                    memory=int(JOB_MEM * 2),
                    runtime=JOB_TIME,
                    cluster_queue=JOB_QUEUE,
                    command_to_run=yld_command,
                )
        else:
            raise ValueError("Must include a YLL version to calculate YLDs")
    if "daly" in stages_to_run:
        if "yld" in stages_to_run:
            if "yll" in stages_to_run:
                daly_args = base_args
                daly_args["save_version"] = daly_save
                daly_args["versions"] = Versions(
                    str(
                        versions.get_version_metadata(base_args["epoch"], "population")
                    ),
                    str(versions.get_version_metadata(base_args["epoch"], "paf")),
                    f"FILEPATH/{yld_save}",
                    f"FILEPATH//{yll_save}",
                )
                daly_args["cur_stage"] = "daly"
                if draws:
                    _ = compute_attributable_burden_code(
                        daly_args, yld_jobs, paf_everywhere=True
                    )
                else:
                    daly_args = args_from_ctx(daly_args)
                    daly_command = (
                        f"{sys.executable} {command} "
                        f"{daly_args} "
                        " compute-attributable-burden "
                    )

                    _ = submit_cluster_job_request(
                        job_name=f"calculate_{desc_string}_daly",
                        cores=JOB_CORES,
                        memory=int(JOB_MEM * 2),
                        runtime=JOB_TIME,
                        cluster_queue=JOB_QUEUE,
                        job_holds=[yld_jid],
                        command_to_run=daly_command,
                    )
            elif "yll" in versions[base_args["epoch"]].keys():
                daly_args = base_args
                daly_args["save_version"] = daly_save
                daly_args["versions"] = Versions(
                    str(
                        versions.get_version_metadata(base_args["epoch"], "population")
                    ),
                    str(versions.get_version_metadata(base_args["epoch"], "paf")),
                    str(versions.get_version_metadata(base_args["epoch"], "yll")),
                    f"FILEPATH/{yld_save}",
                )
                daly_args["cur_stage"] = "daly"
                if draws:
                    _ = compute_attributable_burden_code(
                        death_args, yld_jobs, paf_everywhere=True
                    )
                else:
                    daly_args = args_from_ctx(daly_args)
                    daly_command = (
                        f"{sys.executable} {command} "
                        f"{daly_args} "
                        " compute-attributable-burden "
                    )

                    _ = submit_cluster_job_request(
                        job_name=f"calculate_{desc_string}_daly",
                        cores=JOB_CORES,
                        memory=int(JOB_MEM * 2),
                        runtime=JOB_TIME,
                        cluster_queue=JOB_QUEUE,
                        job_holds=[yld_jid],
                        command_to_run=daly_command,
                    )
            else:
                raise ValueError("Must include a YLL version to calculate DALYs")
        elif "yld" in versions[base_args["epoch"]].keys():
            if "yll" in versions[base_args["epoch"]].keys():
                daly_args = base_args
                daly_args["save_version"] = daly_save
                daly_args["versions"] = Versions(
                    str(
                        versions.get_version_metadata(base_args["epoch"], "population")
                    ),
                    str(versions.get_version_metadata(base_args["epoch"], "paf")),
                    str(versions.get_version_metadata(base_args["epoch"], "yld")),
                    str(versions.get_version_metadata(base_args["epoch"], "yll")),
                )
                daly_args["cur_stage"] = "daly"
                if draws:
                    _ = compute_attributable_burden_code(
                        daly_args, [1], paf_everywhere=True
                    )
                else:
                    daly_args = args_from_ctx(daly_args)
                    daly_command = (
                        f"{sys.executable} {command} "
                        f"{daly_args} "
                        " compute-attributable-burden "
                    )

                    _ = submit_cluster_job_request(
                        job_name=f"calculate_{desc_string}_daly",
                        cores=JOB_CORES,
                        memory=int(JOB_MEM * 2),
                        runtime=JOB_TIME,
                        cluster_queue=JOB_QUEUE,
                        job_holds=[yld_jid],
                        command_to_run=daly_command,
                    )
            else:
                raise ValueError("Must include a YLL version to calculate DALYs")
        else:
            raise ValueError("Must include a YLD version to calculate DALYs")

    logger.info(
        "\n".join(
            [
                "All Jobs Submitted - Save Versions:",
                (f"\tDEATH: {death_save}" if "death" in stages_to_run else ""),
                (f"\tYLLs: {yll_save}" if "yll" in stages_to_run else ""),
                (f"\tYLDs: {yld_save}" if "yld" in stages_to_run else ""),
                (f"\tDALYs: {daly_save}" if "daly" in stages_to_run else ""),
            ]
        )
    )


@cli.command()
@click.pass_context
def compute_attributable_burden(ctx: click.Context) -> None:
    """Compute attributable burden for one stage.

    If --draws command parallelizes by cause.

    Args:
        ctx (click.Context): ctx object passed from cli()
    """
    ctx = ctx_tuple_to_list(ctx)

    if ctx.obj["draws"]:
        compute_attributable_burden_code(ctx.obj, paf_everywhere=True)

    else:
        compute_mean_level_attributable_burden_code(ctx.obj)


@cli.command()
@click.pass_context
def calculate_one_cause_main(ctx: click.Context) -> None:
    """Compute amr attributable burden for one acause.

    Args:
        ctx (click.Context): ctx object passed from cli
    """
    ctx = ctx_tuple_to_list(ctx)

    compute_one_cause(ctx.obj)


if __name__ == "__main__":
    cli()
