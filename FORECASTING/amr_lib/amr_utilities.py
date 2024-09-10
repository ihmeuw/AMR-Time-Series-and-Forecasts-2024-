"""A collection of utility functions for use with AMR forecasting."""

from datetime import datetime
from pathlib import Path
from typing import List, Optional, Union

import click
import xarray as xr
import yaml
from db_queries import get_location_metadata
from fhs_lib_data_transformation.lib.dimension_transformation import expand_dimensions
from fhs_lib_database_interface.lib.constants import AgeConstants, SexConstants
from fhs_lib_file_interface.lib.file_interface import FBDPath
from fhs_lib_file_interface.lib.file_system_manager import FileSystemManager
from fhs_lib_file_interface.lib.os_file_system import OSFileSystem
from fhs_lib_file_interface.lib.version_metadata import VersionMetadata
from fhs_lib_file_interface.lib.versioning import Versions
from fhs_lib_file_interface.lib.xarray_wrapper import open_xr
from fhs_lib_orchestration_interface.lib.slicer import Slicer
from fhs_lib_summary_maker.lib.summary import compute_summary
from fhs_lib_year_range_manager.lib.year_range import YearRange
from tiny_structured_logger.lib import fhs_logging

FileSystemManager.set_file_system(OSFileSystem())

REFERENCE_SCENARIO = 0
COMBINED = 46

NATIONALS_LEVEL = 3
MERGE_COLUMNS = [
    "age_group_id",
    "location_id",
    "sex_id",
    "year_id",
    "scenario",
    "acause",
]
RATE_ONLY = {"age_group_id": [27]}


SCENARIO_FILE_NAMES = {
    -1: "worse",
    0: "reference",
    43: "gram_negative",
    44: "better_care",
    45: "better_care_gram_neg",
    46: "combined",
    52: "paf_carried_forward",
}


TOTAL_DEATH_DIFFERENCE = {
    43: {"total": 0, "reference": 0, "difference": 43},
    45: {"total": 44, "reference": 0, "difference": 43},
}

FHS_LOCATION_SET = 39


def filter_versions(
    versions: dict,
    stage: Optional[List[str]] = None,
    scenario: Optional[List[int]] = None,
    epoch: Optional[List[str]] = None,
    datatype: Optional[List[str]] = None,
) -> dict:
    """Filter versions from versions.yaml file.

    Args:
        versions (dict): versions dictionary returned from get_version_dict
        stage (Optional[List[str]], optional): stage(s) to filter to. Defaults to None.
        scenario (Optional[List[int]], optional): scenario(s) to filter to. Defaults to None.
        epoch (Optional[List[str]], optional): epoch(s) to filter to . Defaults to None.
        datatype (Optional[List[str]], optional): datatype(s) to filter to. Defaults to None.

    Returns:
        dict: Filtered versions dictionary
    """
    if epoch is None:
        epoch = list(versions.keys())
    if scenario is None:
        scenario = list(versions["future" if "future" in epoch else epoch[0]].keys())
    if datatype is None:
        datatype = list(versions[epoch[0]][scenario[0]].keys())
    if stage is None:
        stage = list((versions[epoch[0]][scenario[0]][datatype[0]])[epoch[0]].keys())

    new_dict = {}
    for iter_epoch in epoch:
        new_dict[iter_epoch] = {}
        for iter_scenario in scenario:
            if iter_scenario in versions[iter_epoch]:
                new_dict[iter_epoch][iter_scenario] = {}
                for iter_datatype in datatype:
                    new_dict[iter_epoch][iter_scenario][iter_datatype] = []
                    for iter_stage in stage:
                        new_dict[iter_epoch][iter_scenario][iter_datatype].append(
                            str(
                                flatten_version_dict(
                                    versions,
                                    stage=iter_stage,
                                    scenario=iter_scenario,
                                    epoch=iter_epoch,
                                    datatype=iter_datatype,
                                )["versions"]
                            )
                        )
                    new_dict[iter_epoch][iter_scenario][iter_datatype] = Versions.parse(
                        " ".join(new_dict[iter_epoch][iter_scenario][iter_datatype])
                    )
    return new_dict


def flatten_version_dict(version_dict: dict) -> dict:
    """Turn version dict into a flat dictionary."""
    flat_dict = {}
    if len(version_dict.keys()) > 1:
        raise ValueError
    flat_dict["epoch"] = list(version_dict.keys())[0]

    if len(version_dict[flat_dict["epoch"]]) > 1:
        raise ValueError

    flat_dict["scenario"] = list(version_dict[flat_dict["epoch"]].keys())[0]

    if len(version_dict[flat_dict["epoch"]][flat_dict["scenario"]]) > 1:
        raise ValueError
    flat_dict["datatype"] = list(
        version_dict[flat_dict["epoch"]][flat_dict["scenario"]].keys()
    )[0]

    flat_dict["versions"] = str(
        version_dict[flat_dict["epoch"]][flat_dict["scenario"]][flat_dict["datatype"]]
    )

    return flat_dict


def get_version_dict(yaml_path: str = "versions.py", default_data_dir: int = 7) -> dict:
    """Pull versions from yaml and initialize as Versions.

    Args:
        yaml_path (str, optional): Path to yaml file.
            Must be in directory of amr_utilities.py. Defaults to "versions.py".
        default_data_dir (int, optional): Default data dir for use in Versioning.Versions. Defaults to 7.

    Returns:
        Dict: Dictionary of versions
    """
    cur_dir = Path(__file__).parent.parent.resolve()

    if not (cur_dir / yaml_path).exists():
        if not Path(yaml_path).exists():
            raise FileNotFoundError(
                f"{yaml_path} not found in {str(cur_dir)} or at {str(yaml_path)}"
            )
        else:
            yaml_path = Path(yaml_path)
    else:
        yaml_path = cur_dir / yaml_path

    with open(str(yaml_path)) as yaml_file:
        file_dict = yaml.safe_load(yaml_file)

    for epoch in file_dict:
        for scenario in file_dict[epoch]:
            for datatype in ["attributable", "total", "associated"]:
                if datatype in file_dict[epoch][scenario].keys():
                    file_dict[epoch][scenario][datatype] = Versions.parse(
                        " ".join(file_dict[epoch][scenario][datatype])
                    ).update_default_data_source(default_data_dir)

    return file_dict


def get_one_paf(
    paf_path: Union[FBDPath, Path], acause: str, draw_level: bool = False
) -> xr.DataArray:
    """Calculate PAFs for one cause.

    Args:
        paf_path (Union[FBDPath, Path]): path to PAF version
        acause (str): acause we're pulling
        draw (bool): indicator of if to keep draws

    """
    logger = fhs_logging.get_logger()
    logger.info(f"loading in paf for {acause}")
    acause_da = open_xr(paf_path / f"{acause}.nc")
    if "acause" not in acause_da.coords:
        if "acause" in acause_da.dims:
            acause_da = acause_da.squeeze("acause", drop=True)
        acause_da = expand_dimensions(acause_da, acause=[acause])
    if "draw" not in acause_da.coords:
        if not draw_level:
            logger.warning(
                "No draw coordinate to average across, ignoring pass 'draw' argument"
            )
            draw_level = True
    if not draw_level:
        acause_da = acause_da.chunk(chunks="auto").mean("draw")
    return acause_da


def get_pafs(paf_path: Union[FBDPath, Path], draw_level: bool = False) -> xr.DataArray:
    """Calculate PAFs for all causes.

    Args:
        paf_path (Union[FBDPath, Path]): path to PAF version
        draw (bool): indicator of if to keep draws

    """
    logger = fhs_logging.get_logger()
    all_acause = []
    logger.info("loading in all pafs")
    for acause_path in paf_path.glob("*.nc"):
        all_acause += [get_one_paf(paf_path, acause_path.stem, draw_level)]
    all_acause_da = xr.concat(all_acause, dim="acause")
    return all_acause_da


def get_mean_level_pafs(paf_path: Union[FBDPath, Path]) -> xr.DataArray:
    """Get mean level pafs code for backwards compatibility.

    Args:
        paf_path (Union[FBDPath, Path]): FBDPath to pafs

    Returns:
        xr.DataArray: DataArray with all pafs by acause averaged over draws
    """
    return get_pafs(paf_path, draw_level=False)


def date_now() -> str:
    """Get Date formatted for FHS File Versions.

    Args:
        None

    Returns:
        (str) formatted YYYYMMDD
    """
    return datetime.now().strftime("%Y%m%d")


def filter_dataarrays(
    arg_dict: dict,
    stage: str,
    da: xr.DataArray,
    versions: Versions,
    age_filter: bool = True,
) -> xr.DataArray:
    """Filter passed DataArray based on specifications in arg_dict.

    Args:
        arg_dict (dictionary): arguments to filter xr.DataArray with from cli object
        stage (str): stage filtering data array for (paf, death, yll, etc)
        da (xr.DataArray): DataArray to filter
        versions (Versions): version of data to check for single scenario mode
        age_filter (bool): Flag on whether or not to filter to MD age groups

    Returns:
        xr.DataArray: The passed da with filtered coordinates
    """
    epoch = arg_dict["epoch"]
    if "scenario" in da.coords:
        if versions.get_version_metadata(epoch, stage).scenario is not None:
            da = da.sel(
                scenario=versions.get_version_metadata(epoch, stage).scenario, drop=True
            )
        elif "scenario" in arg_dict:
            da = da.sel(scenario=arg_dict["scenario"], drop=True)
        else:
            da = da.sel(
                scenario=(
                    da.scenario.item()
                    if len(da.scenario.values) == 1
                    else REFERENCE_SCENARIO
                ),
                drop=True,
            )
    if "scenario" in arg_dict:
        da = expand_dimensions(da, scenario=[arg_dict["scenario"]])

    if "statistic" in da.coords:
        da = da.sel(statistic="mean", drop=True)

    if age_filter:
        da = da.sel(age_group_id=AgeConstants.MOST_DETAILED_AGE_GROUP_IDS_PRE_GBD_2020)

    return da


def filter_scenario(
    arg_dict: dict,
    stage: str,
    da: xr.DataArray,
    versions: Optional[Versions] = None,
) -> xr.DataArray:
    """Filter passed DataArray based on specifications in arg_dict.

    Args:
        arg_dict (dictionary): arguments to filter xr.DataArray with from cli object
        stage (str): stage filtering data array for (paf, death, yll, etc)
        da (xr.DataArray): DataArray to filter
        versions (Versions): version of data to check for single scenario mode. Defaults to None.

    Returns:
        xr.DataArray: The passed da with filtered coordinates
    """
    if "scenario" in da.coords:
        if versions is not None:
            epoch = arg_dict["epoch"]
            if versions.get_version_metadata(epoch, stage).scenario is not None:
                da = da.sel(
                    scenario=versions.get_version_metadata(epoch, stage).scenario,
                    drop=True,
                )
        elif "scenario" in arg_dict:
            da = da.sel(scenario=arg_dict["scenario"], drop=True)
        else:
            da = da.sel(
                scenario=(
                    da.scenario.item()
                    if len(da.scenario.values) == 1
                    else REFERENCE_SCENARIO
                ),
                drop=True,
            )

    if "scenario" in arg_dict:
        da = expand_dimensions(da, scenario=[arg_dict["scenario"]])

    return da


def args_from_ctx(arg_object: Union[click.Context, dict]) -> str:
    """Turn click.Context object into a string formatted for command line.

    Args:
        ctx (click.Context): ctx object being passed from cli()

    Returns:
        str formatted '--key value' for single arguments and
        '--key value[0] --key value[1]' for multiple arguments

    """
    if type(arg_object) is not dict:
        arg_object = arg_object.obj
    arg_list = []
    for key, value in arg_object.items():
        if (key == "versions") or isinstance(value, Versions):
            value = tuple([str(x) for x in arg_object["versions"]])
        if key == "save_suffix":
            if value != "":
                arg_list.append(f"--{key.replace('_', '-')} {value}")
        elif key == "rake_subnats":
            arg_list.append(f"--{key.replace('_', '-')} {value}")
        elif isinstance(value, bool):
            if value:
                arg_list.append(f"--{key.replace('_', '-')}")
        elif type(value) in [tuple, list]:
            if len(value) > 0:
                for item in value:
                    arg_list.append(f"--{key.replace('_', '-')} {item}")
        elif isinstance(value, VersionMetadata):
            arg_list.append(f"--{key.replace('_', '-')} {str(value)}")
        elif isinstance(value, YearRange):
            arg_list.append(f"--{key.replace('_', '-')} {str(value)}")
        elif isinstance(value, FBDPath):
            arg_list.append(f"--{key.replace('_', '-')} {str(value.LFN())}")
        elif isinstance(value, Slicer):
            pass
            # arg_list.append(f"--{key.replace('_', '-')} {len(value.draw)}")
        else:
            arg_list.append(f"--{key.replace('_', '-')} {(value)}")
    return " ".join(arg_list)


def ctx_tuple_to_list(ctx: click.Context) -> click.Context:
    """Turn click.Context multiple arguments from tuple to list.

    Args:
        ctx (click.Context): ctx object being passed from cli()

    """
    for key, value in ctx.obj.items():
        if type(value) is tuple:
            ctx.obj[key] = [x for x in value]
    return ctx


def update_totals(
    full_da: xr.DataArray, scenario_run: Optional[int] = None
) -> xr.DataArray:
    """Update total burden with paf level attributable.

    This is for scenarios such as Gram Negative and Worse where
        the scenarios are applied at the PAF level instead
        of through the pipeline so total deaths are never
        updated with the scenarios. This then assumes that
        the total averted deaths for the scenario is
        (SCENARIO ATTRIBUTABLE - REFERENCE ATTRIBUTABLE)
        and updates the total deaths to be
        (TOTAL DEATHS) +
        (SCENARIO ATTRIBUTABLE - REFERENCE ATTRIBUTABLE)

    Args:
        full_da (xr.DataArray): full data returned from load_all_data
        scenario_run (Optional[int]): scenario to update totals for. If none
            does all scenarios in full_da and TOTAL_DEATH_DIFFERENCE
    """
    if scenario_run is not None:
        if scenario_run in TOTAL_DEATH_DIFFERENCE.keys():
            scenario = scenario_run
            total_scenario_da = full_da.sel(
                scenario=scenario,
                datatype="total",
                drop=True,
            )
            reference_attr_da = full_da.sel(
                scenario=TOTAL_DEATH_DIFFERENCE[scenario]["reference"],
                datatype="attributable",
                drop=True,
            )
            difference_attr_da = full_da.sel(
                scenario=TOTAL_DEATH_DIFFERENCE[scenario]["difference"],
                datatype="attributable",
                drop=True,
            )

            # Remove the paf level scenario deaths from total
            total_scenario_da = total_scenario_da + (
                difference_attr_da - reference_attr_da
            )

            # replace just the scenarios totals in full_da
            condition = (full_da.scenario == scenario) & (full_da.datatype == "total")
            full_da = full_da.where(~condition, other=total_scenario_da)
    else:
        for scenario in TOTAL_DEATH_DIFFERENCE:
            if scenario in full_da.scenario.values:
                total_scenario_da = full_da.sel(
                    scenario=scenario,
                    datatype="total",
                    drop=True,
                )
                reference_attr_da = full_da.sel(
                    scenario=TOTAL_DEATH_DIFFERENCE[scenario]["reference"],
                    datatype="attributable",
                    drop=True,
                )
                difference_attr_da = full_da.sel(
                    scenario=TOTAL_DEATH_DIFFERENCE[scenario]["difference"],
                    datatype="attributable",
                    drop=True,
                )

                # Remove the paf level scenario deaths from total
                total_scenario_da = total_scenario_da + (
                    difference_attr_da - reference_attr_da
                )

                # replace just the scenarios totals in full_da
                condition = (full_da.scenario == scenario) & (
                    full_da.datatype == "total"
                )
                full_da = full_da.where(~condition, other=total_scenario_da)

    return full_da


def total_diffs(
    full_da: xr.DataArray, scenario_run: Optional[int] = None
) -> xr.DataArray:
    """Calculate difference of all scenarios compared to reference.

    Args:
        full_da (xr.DataArray): full dataset returned from load_all_data
        scenario_run (Optional[int]): scenario to get the total_diffs for.
            Defaults to  None. If none updates totals for every scenario.
    """
    full_da2 = update_totals(full_da, scenario_run)
    full_da_totals = full_da2

    non_reference_scenario = list(
        set(full_da_totals.scenario.values) - set([REFERENCE_SCENARIO])
    )

    reference_da = full_da_totals.sel(scenario=REFERENCE_SCENARIO, drop=True)

    reference_da = expand_dimensions(reference_da, scenario=non_reference_scenario)

    scenario_da = full_da_totals.sel(scenario=non_reference_scenario)

    return scenario_da - reference_da


def cumulative_diff(
    diff_da: xr.DataArray, sum_range: tuple = (2025, 2050)
) -> xr.DataArray:
    """Sum over passed range and return.

    Args:
        diff_da (xr.DataArray): DataArray with differences
            compared to reference
        sum_range (tuple): Tuple in format (start, end) to sum over
    """
    start, end = sum_range
    cdiff = diff_da.sel(year_id=list(range(start, end + 1))).sum("year_id")
    if "draw" in cdiff.dims:
        cdiff = compute_summary(cdiff.chunk(dict(draw=-1)).chunk(dict(draw=-1)))

    return cdiff


def get_most_detailed(
    da: xr.DataArray,
    release_id: int,
    gbd_2021_ages: bool = False,
    nationals_most_detailed: bool = False,
) -> xr.DataArray:
    """Subset DataArray to it's most detailed demographics.

    Args:
        da (xr.DataArray): DataArray to subset
        release_id (int): Release ID for get_location_metadata
        gbd_2021_ages (bool, optional): Flag to pull GBD 2021 most detailed ages. Defaults to False.
        nationals_most_detailed (bool, optional): Flag to set most detailed locations as nationals. Defaults to False.

    Returns:
        xr.DataArray: DataArray with most detailed demographics
    """
    if gbd_2021_ages:
        ages = AgeConstants.MOST_DETAILED_AGE_GROUP_IDS_GBD_2020_ONWARD
    else:
        ages = AgeConstants.MOST_DETAILED_AGE_GROUP_IDS_PRE_GBD_2020

    locs_metadata = get_location_metadata(
        location_set_id=FHS_LOCATION_SET, release_id=release_id
    )
    if nationals_most_detailed:
        locs = locs_metadata.query("level == @NATIONALS_LEVEL").location_id.to_list()
    else:
        locs = locs_metadata.query("most_detailed == True").location_id.to_list()
    try:
        da = da.sel(age_group_id=ages)
    except KeyError as k:
        missing_ages = list(set(ages) - set(da.age_group_id.values.tolist()))
        print(
            f"Note - subset to md ages failed. Missing age groups: {', '.join([str(x) for x in missing_ages])}"
        )
        da = da.sel(
            age_group_id=list(
                set(ages).intersection(set(da.age_group_id.values.tolist()))
            )
        )
    return da.sel(location_id=locs, sex_id=list(SexConstants.SEX_IDS))
