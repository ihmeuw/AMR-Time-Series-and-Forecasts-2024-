"""Using passed pafs create associated burden forecasts as a ratio of attributable burden.

Example Call:
python make_associated_burden.py \
    --gbd-round-id 7 \
    --release-id 9 \
    --attributable-paf-future VERSION  \
    --attributable-paf-past VERSION \
    --associated-paf-past VERSION \
    --save-version VERSION \
    associated-burden

"""

import sys
from pathlib import Path
from typing import Optional

import click
from amr_lib.amr_utilities import (
    args_from_ctx,
    ctx_tuple_to_list,
    filter_dataarrays,
    filter_scenario,
    get_mean_level_pafs,
    get_most_detailed,
    get_one_paf,
)
from db_queries import get_location_metadata
from fhs_lib_data_transformation.lib.dimension_transformation import expand_dimensions
from fhs_lib_file_interface.lib.file_interface import FBDPath
from fhs_lib_file_interface.lib.file_system_manager import FileSystemManager
from fhs_lib_file_interface.lib.os_file_system import OSFileSystem
from fhs_lib_file_interface.lib.version_metadata import VersionMetadata
from fhs_lib_file_interface.lib.versioning import Versions
from fhs_lib_file_interface.lib.xarray_wrapper import save_xr
from fhs_lib_orchestration_interface.lib.submit_cluster_job import (
    submit_cluster_job_request,
)
from fhs_lib_year_range_manager.lib.year_range import YearRange
from tiny_structured_logger.lib import fhs_logging

FileSystemManager.set_file_system(OSFileSystem())

FHS_LOCATION_SET = 39
FHS_CAUSE_SET = 6
EPOCHS = ["past", "future"]
SCENARIO_NAMES = {
    -1: "Pessimistic scenario (15th percentile)",
    0: "Reference",
    43: "Gram Negative",
    44: "Better Access to Antibiotics",
    46: "Combined (Gram Negative, Better Access, WaSH Elimination, Vaccine Scale-Up)",
}


SCENARIO_FILE_NAMES = {
    -1: "worse",
    0: "reference",
    43: "gram_negative",
    44: "better_access",
    46: "combined",
}

SEX_NAMES = {1: "Male", 2: "Female", 3: "Both"}

TOTAL_DEATH_DIFFERENCE = {
    -1: {"total": 0, "reference": 0, "difference": -1},
    43: {"total": 0, "reference": 0, "difference": 43},
    46: {"total": 46, "reference": 46, "difference": 43},
}

YEARS = YearRange.parse_year_range("1990:2022:2050")

# This is the query to get all of the acause's we'll be reporting
# to the AMR team. Note _oth_pand is not forecast so will not be
# included.
CODEBOOK_CAUSES = "(level <= 2) and acause != '_oth_pand'"
REPORTING_LOCATIONS = "level <= 1"

REFERENCE_SCENARIO = 0
GRAM_NEGATIVE = 43
COMBINED = 46
WORSE = -1
BETTER = 1


def make_associated_burden_one_cause(arg_dict: dict) -> None:
    """Make Associated burden for one acause at the draw level.

    Args:
        arg_dict (dict): Dictionary of args from CLI
    """
    logger = fhs_logging.get_logger()
    gbd_round_id = arg_dict["gbd_round_id"]
    release_id = arg_dict["release_id"]
    attributable_paf_future = arg_dict["attributable_paf_future"]
    attributable_paf_past = arg_dict["attributable_paf_past"]
    associated_paf_past = arg_dict["associated_paf_past"]
    save_version = arg_dict["save_version"]
    scenario = arg_dict["scenario"]
    acause = arg_dict["acause"]

    logger.info("Loading in all past files and filtering to national level")
    attributable_paf_future_da = get_one_paf(
        attributable_paf_future.data_dir(), acause=acause, draw_level=True
    )
    attributable_paf_past_da = get_one_paf(
        attributable_paf_past.data_dir(), acause=acause, draw_level=True
    )
    associated_paf_past_da = get_one_paf(
        associated_paf_past.data_dir(), acause=acause, draw_level=True
    )

    attributable_paf_future_da = get_most_detailed(
        attributable_paf_future_da, release_id, nationals_most_detailed=True
    )
    attributable_paf_past_da = get_most_detailed(
        attributable_paf_past_da, release_id, nationals_most_detailed=True
    )
    associated_paf_past_da = get_most_detailed(
        associated_paf_past_da, release_id, nationals_most_detailed=True
    )

    attributable_paf_future_da = filter_scenario(
        arg_dict, "paf", attributable_paf_future_da
    )
    attributable_paf_past_da = filter_scenario(
        arg_dict, "paf", attributable_paf_past_da
    )
    associated_paf_past_da = filter_scenario(arg_dict, "paf", associated_paf_past_da)

    associated_ratio = associated_paf_past_da.sel(
        year_id=YEARS.past_end, drop=True
    ) / attributable_paf_past_da.sel(year_id=YEARS.past_end, drop=True)

    associated_paf_future = attributable_paf_future_da * associated_ratio

    save_xr(
        associated_paf_future,
        save_version / f"{acause}.nc",
        metric="percent",
        space="identity",
    )


def make_mean_level_associated_burden(arg_dict: dict) -> None:
    """Make associated burden forecasts using passed arguments.

    Args:
        arg_dict (dict): Objects passed to cli.
    """
    logger = fhs_logging.get_logger()
    gbd_round_id = arg_dict["gbd_round_id"]
    release_id = arg_dict["release_id"]
    attributable_paf_future = arg_dict["attributable_paf_future"]
    attributable_paf_past = arg_dict["attributable_paf_past"]
    associated_paf_past = arg_dict["associated_paf_past"]
    save_version = arg_dict["save_version"]
    scenario = arg_dict["scenario"]

    location_metadata = get_location_metadata(
        release_id=release_id, location_set_id=FHS_LOCATION_SET
    )

    md_locs = location_metadata.query("level == 3").location_id.to_list()

    logger.info("Loading in all past files and filtering to national level")

    attributable_paf_future_da = get_mean_level_pafs(
        attributable_paf_future.data_dir()
    ).sel(location_id=md_locs, sex_id=[1, 2])
    attributable_paf_past_da = get_mean_level_pafs(
        attributable_paf_past.data_dir()
    ).sel(location_id=md_locs)
    associated_paf_past_da = get_mean_level_pafs(associated_paf_past.data_dir()).sel(
        location_id=md_locs, sex_id=[1, 2]
    )

    filter_dict = arg_dict

    filter_dict["epoch"] = "future"
    filter_versions = Versions(str(attributable_paf_future))
    attributable_paf_future_da = filter_dataarrays(
        filter_dict, "paf", attributable_paf_future_da, filter_versions, age_filter=True
    )

    filter_dict["epoch"] = "past"

    filter_versions = Versions(str(attributable_paf_past))
    attributable_paf_past_da = filter_dataarrays(
        filter_dict, "paf", attributable_paf_past_da, filter_versions, age_filter=True
    )

    filter_versions = Versions(str(associated_paf_past))
    associated_paf_past_da = filter_dataarrays(
        filter_dict, "paf", associated_paf_past_da, filter_versions, age_filter=True
    )

    associated_ratio = associated_paf_past_da.sel(
        year_id=YEARS.past_end, drop=True
    ) / attributable_paf_past_da.sel(year_id=YEARS.past_end, drop=True)

    if "scenario" in associated_ratio.coords:
        associated_ratio = associated_ratio.sel(
            scenario=associated_ratio.scenario.item(), drop=True
        )
    associated_paf_future = attributable_paf_future_da * associated_ratio
    associated_paf_future = expand_dimensions(
        associated_paf_future, scenario=[scenario]
    )
    for acause in associated_paf_future.acause.values:
        acause_da = associated_paf_future.sel(acause=acause, drop=True)
        save_xr(
            acause_da, save_version / f"{acause}.nc", metric="percent", space="identity"
        )


with_release_id = click.option(
    "--release-id",
    type=int,
    required=False,
    default=9,
    help="Release ID for pulling location sets",
)

with_gbd_round_id = click.option(
    "--gbd-round-id", type=int, default=7, help="GBD Round ID"
)


with_attributable_burden = click.option(
    "--attributable_burden",
    "-a",
    type=str,
    multiple=True,
    callback=Versions.parse_cli_versions_callback,
    help="Attributable burden versions to use",
)
with_attributable_paf_future = click.option(
    "--attributable-paf-future",
    type=str,
    required=True,
    help="Future paf to apply ratio to",
)

with_attributable_paf_past = click.option(
    "--attributable-paf-past",
    type=str,
    required=True,
    help="past paf to calculate ratio with",
)
with_associated_paf_past = click.option(
    "--associated-paf-past",
    type=str,
    required=True,
    help="Past paf to calculate ratio with",
)

with_scenario = click.option(
    "--scenario",
    "-s",
    type=int,
    default=REFERENCE_SCENARIO,
    help="scenario ID to pull from PAFs",
)

with_save_version = click.option(
    "--save-version", type=str, help="Version to save file as"
)

with_acause = click.option("--acause", type=str, required=False, help="Acause to pull")


@click.group()
@with_release_id
@with_gbd_round_id
@with_associated_paf_past
@with_save_version
@with_attributable_paf_future
@with_attributable_paf_past
@with_scenario
@with_acause
@click.pass_context
def cli(
    ctx: click.Context,
    gbd_round_id: int,
    release_id: int,
    attributable_paf_future: str,
    attributable_paf_past: str,
    associated_paf_past: str,
    scenario: str,
    acause: Optional[str],
    save_version: str,
) -> None:
    """Initialize cli function.

    Args:
        ctx (click.Context): ctx object passed from Click
        gbd_round_id (int): GBD Round ID for use in file location
        release_id (int): Release ID for us in db_queries functions
        attributable_paf_future (str): Forecasted PAFs of Attributable Burden
        attributable_paf_past (str): ETL'd PAFs of Attributable Burden
        associated_paf_past (str): ETL'd PAFs of Associated Burden
        scenario (str): Scenario ID of PAFs and to Save version as
        acause (Optional[str]): Acause to pull for parallelization
        save_version (str): Version to Save Forecasted Associated Burden PAFs under
    """
    ctx.obj = {
        "gbd_round_id": gbd_round_id,
        "release_id": release_id,
        "attributable_paf_future": VersionMetadata.parse_version(
            f"FILEPATH/{attributable_paf_future}"
        ),
        "attributable_paf_past": VersionMetadata.parse_version(
            f"FILEPATH/{attributable_paf_past}"
        ),
        "associated_paf_past": VersionMetadata.parse_version(
            f"FILEPATH/{associated_paf_past}"
        ),
        "scenario": scenario,
        "acause": acause,
        "save_version": FBDPath(f"FILEPATH/{save_version}"),
    }


@cli.command()
@click.pass_context
def mean_level_associated_burden(ctx: click.Context) -> None:
    """Make mean level future associated burden based on past files.

    Args:
        ctx (click.Context): ctx object passed from cli()
    """
    ctx = ctx_tuple_to_list(ctx)

    make_mean_level_associated_burden(ctx.obj)


@cli.command()
@click.pass_context
def draw_level_associated_burden(ctx: click.Context) -> None:
    """Make draw level future associated burden based on past files in parallel by acause.

    Args:
        ctx (click.Context): ctx object passed from cli()
    """
    ctx = ctx_tuple_to_list(ctx)
    paf_path = ctx.obj["attributable_paf_future"]

    ctx.obj["attributable_paf_future"] = ctx.obj["attributable_paf_future"].version
    ctx.obj["attributable_paf_past"] = ctx.obj["attributable_paf_past"].version
    ctx.obj["associated_paf_past"] = ctx.obj["associated_paf_past"].version
    ctx.obj["save_version"] = ctx.obj["save_version"].stem

    python_call = sys.executable
    filepath = str(Path(__file__).absolute())

    args = args_from_ctx(ctx.obj)

    for acause_path in paf_path.data_dir().glob("*.nc"):
        function_call = f"""
        {python_call} {filepath} {args} --acause {acause_path.stem} one-cause-main
        """
        submit_cluster_job_request(
            f"{acause_path.stem}_associated_burden",
            memory=50,
            cores=5,
            command_to_run=function_call,
            cluster_queue="all.q",
        )


@cli.command()
@click.pass_context
def one_cause_main(ctx: click.Context) -> None:
    """Make draw level future associated burden for one cause.

    Args:
        ctx (click.Context): ctx object passed from cli()
    """
    ctx = ctx_tuple_to_list(ctx)

    make_associated_burden_one_cause(ctx.obj)


if __name__ == "__main__":
    cli()
