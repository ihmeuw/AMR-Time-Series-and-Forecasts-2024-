"""This script uses a cubic spline to interpolate syndrome CFRs by HAQi (by age).

The interpolated values are used to calculate reduction scalars used in `calculate_better_access_adjusted_mx.py`.

"""

import click
import xarray as xr
from db_queries import get_location_metadata
from fhs_lib_file_interface.lib.file_interface import FBDPath
from fhs_lib_file_interface.lib.file_system_manager import FileSystemManager
from fhs_lib_file_interface.lib.os_file_system import OSFileSystem
from fhs_lib_file_interface.lib.xarray_wrapper import open_xr, save_xr
from fhs_lib_year_range_manager.lib.year_range import YearRange

FileSystemManager.set_file_system(OSFileSystem())


FHS_LOCATION_SET_ID = 39


@click.command()
@click.option("--input-version", type=str, help="Version to pull CFRs from")
@click.option("--output-version", type=str, help="Version to Save Interpolated CFRs to")
@click.option("--haq-version", type=str, help="HAQ version to pull from")
@click.option("--gbd-round-id", type=int, help="gbd round id for saving")
@click.option("--release-id", type=int, help="release id for db_queries")
@click.option("--years", type=YearRange.parse_year_range, help="Years")
@click.option("--target-year", type=int, help="Target year for interpolation")
@click.option("--percentile", type=float, default=0.85, help="Percentile for goal HAQ")
def interpolate(
    input_version: str,
    output_version: str,
    haq_version: str,
    gbd_round_id: int,
    release_id: int,
    years: YearRange,
    target_year: int,
    percentile: float,
) -> None:
    """Interpolate ETL'd CFR by HAQi.

    Args:
        input_version (str): CFRs ETL'd from AMR Team
        output_version (str): Location to Save interpolated CFRs
        haq_version (str): HAQ version to use
        gbd_round_id (int): GBD Round ID to save to
        release_id (int): Release ID for db_queries
        years (YearRange): Years object
        target_year (int): Target year to interpolate to
        percentile (float): Percentile to interpolate to
    """
    cfr_path = FBDPath(f"FILEPATH/{input_version}")
    output_path = FBDPath(f"FILEPATH/{output_version}")
    haq_path = FBDPath("FILEPATH/haq.nc")

    locs = get_location_metadata(
        release_id=release_id, location_set_id=FHS_LOCATION_SET_ID
    ).query("level==3")

    percentile_year = years.past_end
    past_end = years.past_end
    past_years = years.past_years

    haq = (
        open_xr(haq_path)
        .sel(statistic="mean", location_id=locs.location_id.values)
        .drop(["statistic", "age_group_id", "sex_id"])
        .squeeze(["sex_id", "age_group_id"])
        .rename("haq")
    )

    cfr_files = [fpath for fpath in cfr_path.glob("*.nc")]

    # Calculate reduction scalars
    syndrome_cfr_by_haq = []
    syndrome_relative_reduction_scalars = []
    for i, file in enumerate(cfr_files):

        print(file.stem)

        cfr = open_xr(file).sel(year_id=past_years).rename("cfr")

        if i == 0:
            haq = haq.sel(location_id=cfr.location_id.values)
            haq_df = haq.to_dataframe().reset_index()
            eighty_fifth = haq.sel(year_id=percentile_year).quantile(percentile).item()

        cfr_df = cfr.to_dataframe().reset_index()
        both_df = haq_df.merge(cfr_df)

        cfr_by_haq = (
            both_df.drop(columns=["year_id", "location_id"])
            .set_index(["haq", "age_group_id", "syndrome"])
            .to_xarray()["cfr"]
        )

        cfr_interpolated = cfr_by_haq.interp(
            haq=eighty_fifth, method="cubic"
        ).expand_dims("haq")

        cfr_by_haq_with_target = xr.concat([cfr_by_haq, cfr_interpolated], dim="haq")
        syndrome_cfr_by_haq.append(cfr_by_haq_with_target)

        cfr_past_end = cfr.sel(year_id=past_end, drop=True)

        reduction = (
            (cfr_interpolated / cfr_past_end)
            .drop("haq")
            .squeeze("haq")
            .assign_coords(year_id=target_year)
            .expand_dims("year_id")
        )

        reduction_scalar = 1 - reduction.clip(max=1, min=0)

        syndrome_relative_reduction_scalars.append(reduction_scalar)

    # Save reduction scalars
    for da in syndrome_relative_reduction_scalars:
        syndrome = da.syndrome.item()
        save_xr(da, output_path / f"{syndrome}.nc", metric="rate", space="identity")
        print(f"Saved {syndrome} to {output_path}")


if __name__ == "__main__":
    interpolate()
