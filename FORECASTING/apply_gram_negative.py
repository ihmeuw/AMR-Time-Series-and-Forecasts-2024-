"""Apply Gram negative fractions to forecasted PAFs.

This script applies gram negative fractions to forecasted PAFS in a few steps:
1. Read in past gram negative fractions
2. Read in future PAF REFs
3. Multiply 2021 past gram negative by future PAF REF
4. Multiply 1 - (2021 past gram negative) by future PAF REF
5. Create Gram negative scenario where 2036 is 50% of 2021 value
6. Add back on Gram Positive

The files are saved in the following structure:
- gram_negative_ref: REFs with gram negative fractions
- gram_positive: REFs with gram positive fractions
- gram_negative_scen: Scenario with gram negative fractions
- final: Scenario with gram positive and gram negative fractions

Final is saved in the main version and the other three are saved in subfolders.

Checks you should run:
1. 2036 value in final_future should be half that of 2021
2. 2037 on in final_future is constant
3. 2021 in new_scenario is equal to original future
4. Other values are different
5. Constant decrease in final_future

Example usage:

bash::

python apply_gram_negative.py \
    apply-gram-negative \
    --gbd_round_id 6 \
    --gram-negative-path FILEPATH \
    --years 1990:2021:2050 \
    --paf-version FILEPATH \
    --output-version FILEPATH \
    --target-year 2036 \
    --target-prop 0.5
"""

import click
import xarray as xr
from fhs_lib_data_transformation.lib.dimension_transformation import expand_dimensions
from fhs_lib_database_interface.lib.query.location import get_location_set
from fhs_lib_file_interface.lib.file_interface import FBDPath
from fhs_lib_file_interface.lib.file_system_manager import FileSystemManager
from fhs_lib_file_interface.lib.os_file_system import OSFileSystem
from fhs_lib_file_interface.lib.version_metadata import VersionMetadata
from fhs_lib_file_interface.lib.xarray_wrapper import open_xr, save_xr
from fhs_lib_year_range_manager.lib.year_range import YearRange

LOCATION_QUERY = "level == 3"
REFERENCE = 0
GRAM_NEGATIVE = 43

FileSystemManager.set_file_system(OSFileSystem())


@click.command()
@click.option(
    "--gbd-round-id",
    type=int,
    default=6,
    help="gbd_round ID to pull fhs location set with",
)
@click.option(
    "--gram-negative-path",
    type=FBDPath,
    help="path to ETL'd gram negative fractions",
)
@click.option(
    "--years", type=YearRange.parse_year_range, help="Years to use for calculations"
)
@click.option(
    "--paf-version",
    type=VersionMetadata.default_parser(default_epoch="future", default_stage="paf"),
    help="PAF put through GenEM version to use as input",
)
@click.option(
    "--output-version",
    type=VersionMetadata.default_parser(default_epoch="future", default_stage="paf"),
    help="Location to save as output",
)
@click.option(
    "--target-year",
    type=int,
    default=2036,
    help="target year to reach target_prop of 2021 value",
)
@click.option(
    "--target-prop",
    type=float,
    default=0.5,
    help="Target proportion of 2021 value to reach by target_year",
)
def apply_gram_negative(
    gbd_round_id: int,
    gram_negative_path: FBDPath,
    years: YearRange,
    paf_version: VersionMetadata,
    output_version: VersionMetadata,
    target_year: int,
    target_prop: float,
) -> None:
    """Apply gram negative fraction to passed paf_version.

    Args:
        gbd_round_id (int): gbd_round ID to pull fhs location set with
        gram_negative_path (FBDPath): path to ETL'd gram negative fractions
        years (YearRange): Years to use for calculations
        paf_version (VersionMetadata): PAF put through GenEM version to use as input
        output_version (VersionMetadata): Location to save as output
        target_year (int): target year to reach target_prop of 2021 value
        target_prop (float): Target proportion of 2021 value to reach by target_year
    """
    paf_version = paf_version.default_data_source(gbd_round_id)
    output_version = output_version.default_data_source(gbd_round_id)
    loc_6_full = get_location_set(gbd_round_id=gbd_round_id)
    loc_countries = loc_6_full.query(LOCATION_QUERY).location_id.values
    fractions = open_xr(gram_negative_path)
    fractions = fractions.sel(
        year_id=years.past_end, location_id=loc_countries, drop=True
    )
    for filepath in paf_version.data_dir().glob("*.nc"):
        entity = filepath.stem
        ref = open_xr(filepath)
        attrs = ref.attrs
        ref = ref.sel(scenario=REFERENCE, drop=True)
        gram_neg_ref = ref * fractions
        gram_pos = ref * (1 - fractions)

        save_xr(
            gram_neg_ref,
            output_version.data_dir() / f"gram_negative_ref/{entity}.nc",
            **attrs,
        )
        save_xr(
            gram_pos, output_version.data_dir() / f"gram_positive/{entity}.nc", **attrs
        )

        past_last_year = gram_neg_ref.sel(year_id=years.past_end, drop=True)
        past_with_future_data = [gram_neg_ref.sel(year_id=[years.past_end])]
        scenario_year_range = target_year - years.past_end
        end_data = past_last_year * target_prop
        rate_change = (past_last_year - end_data) / scenario_year_range

        for ix, year in enumerate(range(years.forecast_start, target_year + 1)):
            subtracted_values = (ix + 1) * rate_change
            newest_year_value = past_last_year - subtracted_values
            past_with_future = expand_dimensions(newest_year_value, year_id=[year])
            past_with_future_data.append(past_with_future)

        past_with_future_da = xr.concat(past_with_future_data, dim="year_id")

        # Take the values from 2036 and hold them constant to 2050
        past_with_future_end_year = past_with_future_da.sel(year_id=target_year).drop(
            "year_id"
        )
        past_with_future_constant = expand_dimensions(
            past_with_future_end_year,
            year_id=range(target_year + 1, years.forecast_end + 1),
        )
        final_future = xr.concat(
            [past_with_future_da, past_with_future_constant], dim="year_id"
        )
        final_future = expand_dimensions(final_future, scenario=[GRAM_NEGATIVE])
        save_xr(
            final_future,
            output_version.data_path() / f"gram_negative_scen/{entity}.nc",
            **attrs,
        )

        new_scenario = final_future + gram_pos

        save_xr(new_scenario, output_version.data_path() / f"{entity}.nc", **attrs)


if __name__ == "__main__":
    apply_gram_negative()
