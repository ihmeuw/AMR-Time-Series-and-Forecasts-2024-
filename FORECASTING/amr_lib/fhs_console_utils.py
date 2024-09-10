"""Helper functions to run aggregator and summary-maker with submit_cluster_job_request."""

import sys
from pathlib import Path
from typing import List, Optional

from amr_lib.amr_utilities import args_from_ctx
from fhs_lib_file_interface.lib.versioning import Versions
from fhs_lib_orchestration_interface.lib.submit_cluster_job import (
    submit_cluster_job_request,
)
from fhs_lib_year_range_manager.lib.year_range import YearRange

POPULATION = "population"


def run_aggregator_console(
    stage: str,
    gbd_round_id: int,
    years: str,
    draws: int,
    epoch: str,
    output_version: Optional[str] = None,
    output_scenario: Optional[int] = None,
    version_to_aggregate: Optional[str] = None,
    population_version: Optional[str] = None,
    versions: Optional[Versions] = None,
    aggregation_round_id: Optional[int] = None,
    holds: List = [1],
    extra_flags: List = [],
    memory: int = 100,
    queue: str = "all",
    scenario: int = None,
    runtime: str = "04:00:00",
    job_name: str = "aggregate",
    log_level: str = "INFO",
    provenance: bool = False,
    save_number: bool = True,
    save_percent: bool = False,
    cause_aggregation: bool = True,
    demographic_aggregation: bool = True,
    rake_subnats: bool = False,
    national_only: bool = True,
    custom_cause_hierarchies: List = [],
    custom_location_set_ids: List = [],
) -> int:
    """Run fhs_lib_data_aggregation_console aggregate-whole-stage.

    Args:
        stage (str): stage to aggregate
        gbd_round_id (int): gbd_round_id for FBDpath and used for
            hierarchies if aggregation_round_id is not passed
        years (str): years for passing to aggregator argument.
            Should be in format "YYYY:YYYY:YYYY"
        draws (int): draws for passing to aggregator. Should
            be number of draws in files.
        epoch (str): Epoch we are aggregating. Can be past or future.
        version_to_aggregate (Optional[str]): Optional version to aggregate.
            Not needed if passing in a Versions object
        population_version (Optional[str]): Optional population to use for
            aggregates. Not needed if passing in a Versions object
        versions (Optional[Versions]): Optional Versions to use for aggregating.
            Not needed if passing in version_to_aggregate and population_version
        aggregation_round_id (Optional[int]): Optional gbd_round_id to use to
            get hierarchies and metadata for aggregation. For use with round
            shifting.
        holds (List, optional): Job dependencies to pass into
            submit_cluster_job_request. Defaults to [1].
        extra_flags (List, optional): extra_flags to pass into
            submit_cluster_job_request. Defaults to [].
        memory (int, optional): Memory (in GB) of aggregation Job. Defaults to 100.
        queue (str, optional): Queue to run aggregation job on. Defaults to "all".
        scenario (int, optional): Scenario to pass into aggregator. Defaults to 0.
        runtime (str, optional): Max runtime to allow aggregation job. Defaults to "04:00:00".
        job_name (str, optional): Job name to pass into
            submit_cluster_job_request. Defaults to "aggregate".
        log_level (str, optional): Log level to pass to aggregator. Defaults to "INFO".
        provenance (bool, optional): Flag to run aggregator with provenance. Defaults to False.
        save_number (bool, optional): Flag to have aggregator save number as well as rate. Defaults to True.
        save_percent (bool, optional): Flag to have aggregator save percent. Defaults to False.
        cause_aggregation (bool, optional): Flag to have aggregator aggregate causes. Defaults to True.
        demographic_aggregation (bool, optional): Flag to have aggregator aggregate demographics. Defaults to True.
        rake_subnats (bool, optional): Flag to have aggregator rake subnats. Defaults to False.
        national_only (bool, optional): Flag to have aggregator ignore subnats. Defaults to True.
        custom_cause_hierarchies (List, optional): List of custom cause hierarchy IDs. Defaults to [].
        custom_location_set_ids (List, optional): List of custom location set IDs. Defaults to [].

    Raises:
        ValueError: if versions isn't passed and missing version_to_aggregate and/or
            population_version

    Returns:
        int: Job ID returned from submit_cluster_job_request
    """
    if versions is None:
        if (version_to_aggregate is None) or (population_version is None):
            raise ValueError(
                "If versions object is not past, must specify version_to_aggregate and population_version"
            )

    if aggregation_round_id is None:
        aggregation_round_id = gbd_round_id

    agg_args = {
        "gbd_round_id": aggregation_round_id,
        "stage": stage,
        "years": years,
        "draws": draws,
        "log_level": log_level,
        "provenance": provenance,
        "save_number": save_number,
        "save_percent": save_percent,
        "cause_aggregation": cause_aggregation,
        "demographic_aggregation": demographic_aggregation,
        "rake_subnats": rake_subnats,
        "national_only": national_only,
        "custom_cause_hierarchies": custom_cause_hierarchies,
        "custom_location_set_ids": custom_location_set_ids,
        "cluster_queue": f"{queue}.q",
    }
    if output_version is not None:
        agg_args["output_version"] = output_version
    else:
        agg_args["output_version"] = (
            f"FILEPATH/{version_to_aggregate}_agg"
            if versions is None
            else versions.get_version_metadata(epoch, stage).with_version(
                f"{versions.get_version_metadata(epoch, stage).version}_agg"
            )
        )
    if output_scenario is not None:
        agg_args["output_scenario"] = output_scenario
    agg_args["versions"] = (
        versions
        if versions is not None
        else Versions(
            str(f"FILEPATH/{version_to_aggregate}"),
            str(f"FILEPATH/{population_version}"),
        )
    )

    agg_args = args_from_ctx(agg_args)

    executable_path = str(
        Path(sys.executable).parent / "fhs_lib_data_aggregation_console"
    )
    aggregation_call = f"{executable_path} aggregate-whole-stage {agg_args}"

    return submit_cluster_job_request(
        job_name=job_name,
        memory=memory,
        cores=1,
        command_to_run=aggregation_call,
        runtime=runtime,
        cluster_queue=f"{queue}.q",
        job_holds=holds,
        extra_flags=extra_flags,
    )


def run_summary_maker_console(
    stage: str,
    gbd_round_id: int,
    years: str,
    epoch: str,
    version_to_summarize: Optional[str] = None,
    versions: Optional[Versions] = None,
    holds: List = [1],
    extra_flags: List = [],
    memory: int = 100,
    queue: str = "all",
    runtime: str = "04:00:00",
    job_name: str = "summarize",
    log_level: str = "INFO",
    national_only: bool = False,
    combine_past_and_future: bool = False,
    compute_demographic_aggregates: bool = False,
    entity_type: str = "acause",
    input_entities_aggregated: str = True,
    no_missing_entity_check: bool = True,
) -> int:
    """Run fhs_lib_summary_maker_console parallelize-entities with specified args.

    Args:
        stage (str): stage to aggregate
        gbd_round_id (int): gbd_round_id for FBDpath and used for
            hierarchies if aggregation_round_id is not passed
        years (str): years for passing to aggregator argument.
            Should be in format "YYYY:YYYY:YYYY"
        epoch (str): Epoch we are aggregating. Can be past or future.
        version_to_summarize (Optional[str]): Optional version to make summary of.
            Not needed if passing in a Versions object.
        versions (Optional[Versions]): Optional Versions to use with summary maker.
            Not needed if passing in version_to_summarize.
        holds (List, optional): Job dependencies to pass into
            submit_cluster_job_request. Defaults to [1].
        extra_flags (List, optional): extra_flags to pass into
            submit_cluster_job_request. Defaults to [].
        memory (int, optional): Memory (in GB) of aggregation Job. Defaults to 100.
        queue (str, optional): Queue to run aggregation job on. Defaults to "all".
        runtime (str, optional): Max runtime to allow aggregation job. Defaults to "04:00:00".
        job_name (str, optional): Job name to pass into
            submit_cluster_job_request. Defaults to "aggregate".
        log_level (str, optional): Log level to pass to aggregator. Defaults to "INFO".
        national_only (bool, optional): Flag to exclude subnationals. Defaults to False.
        combine_past_and_future (bool, optional): Flag to combine past and future. Defaults to False.
        compute_demographic_aggregates (bool, optional): Flag to aggregate locations. Defaults to False.
        entity_type (str, optional): Entity type of stage. Defaults to "acause".
        input_entities_aggregated (str, optional): Flag indicating the the input is aggregated. Defaults to True.
        no_missing_entity_check (bool, optional): Flag to not fail if any most detailed
            entities are missing. Defaults to True.

    Raises:
        ValueError: if versions isn't passed and missing version_to_summarize is missing

    Returns:
        int: Job ID returned from submit_cluster_job_request
    """
    if versions is None:
        if version_to_summarize is None:
            raise ValueError(
                "If versions object is not past, must specify version_to_summarize"
            )

    summary_args = {
        "gbd_round_id": gbd_round_id,
        "years": years,
        "stage": stage,
        "log_level": log_level,
        "national_only": national_only,
        "combine_past_and_future": combine_past_and_future,
        "compute_demographic_aggregates": compute_demographic_aggregates,
        "entity_type": entity_type,
        "input_entities_aggregated": input_entities_aggregated,
        "no_missing_entity_check": no_missing_entity_check,
        "cluster_queue": f"{queue}.q",
    }

    summary_args["versions"] = (
        versions
        if versions is not None
        else Versions(str(f"7/{epoch}/{stage}/{version_to_summarize}"))
    )

    command_args = args_from_ctx(summary_args)

    executable_path = str(Path(sys.executable).parent / "fhs_lib_summary_maker_console")
    summarize_call = f"{executable_path} parallelize-entities {command_args}"

    return submit_cluster_job_request(
        job_name=job_name,
        memory=memory,
        cores=1,
        command_to_run=summarize_call,
        runtime=runtime,
        cluster_queue=f"{queue}.q",
        job_holds=holds,
        extra_flags=extra_flags,
    )


def run_yll_console(
    versions: Versions,
    output_scenario: int,
    gbd_round_id: int,
    draws: int,
    reference_lex_gbd_round_id: int,
    years: YearRange,
    save_number: bool = True,
    log_level: str = "INFO",
    provenance: bool = False,
    holds: List = [1],
    extra_flags: List = [],
    memory: int = 100,
    queue: str = "all",
    runtime: str = "04:00:00",
    job_name: str = "calculate_ylls",
) -> int:
    """Run fhs_pipeline_yll_console parallelize with specified arguments.

    Args:
        versions (Versions): Versions to pass into co
        output_scenario (int): Output scenario for YLLs
        gbd_round_id (int): GBD Round ID for aggregation
        draws (int): Draws in Input and output files
        reference_lex_gbd_round_id (int): GBD Round ID of reference LEX
        years (YearRange): Years to pass into command
        save_number (bool, optional): Flag to save numbers. Defaults to True.
        log_level (str, optional): Log level for pipeline. Defaults to "INFO".
        provenance (bool, optional): Flag to run with provenance. Defaults to False.
        holds (List, optional): Optional Holds to pass to submit_cluster_job_request. Defaults to [1].
        extra_flags (List, optional): Custom arguments to pass to submit_cluster_job_request. Defaults to [].
        memory (int, optional): Memory to pass to submit_cluster_job_request. Defaults to 100.
        queue (str, optional): Queue to pass to submit_cluster_job_request. Defaults to "all".
        runtime (str, optional): Runtime to pass to submit_cluster_job_request. Defaults to "04:00:00".
        job_name (str, optional): Job name to pass to submit_cluster_job_request. Defaults to "calculate_ylls".

    Returns:
        int: Job ID for submitted request.
    """
    yll_args = {
        "log_level": log_level,
        "gbd_round_id": 6,  # gbd_round_id,
        "draws": draws,
        "save_number": save_number,
        "years": years,
        "output_scenario": output_scenario,
        "reference_lex_gbd_round_id": reference_lex_gbd_round_id,
        "provenance": provenance,
        "versions": versions,
    }

    command_args = args_from_ctx(yll_args)

    executable_path = str(Path(sys.executable).parent / "fhs_pipeline_yll_console")
    summarize_call = f"{executable_path} parallelize {command_args}"

    return submit_cluster_job_request(
        job_name=job_name,
        memory=memory,
        cores=1,
        command_to_run=summarize_call,
        runtime=runtime,
        cluster_queue=f"{queue}.q",
        job_holds=holds,
        extra_flags=extra_flags,
    )


def run_daly_console(
    gbd_round_id: int,
    years: YearRange,
    draws: int,
    versions: Versions,
    output_scenario: Optional[int] = None,
    inputs_aggregated: bool = False,
    save_number: bool = True,
    provenance: bool = False,
    holds: List = [1],
    extra_flags: List = [],
    memory: int = 100,
    queue: str = "all",
    runtime: str = "04:00:00",
    job_name: str = "calculate_dalys",
) -> int:
    """Run fhs_pipeline_dalys_console parallelize-by-cause with specified arguments.

    Args:
        gbd_round_id (int): GBD Round ID to pass into command
        years (YearRange): Years to pass into command
        draws (int): Draws of input and output files
        versions (Versions): Versions - must include YLL, YLDs, DALY and population
        output_scenario (Optional[int]): Scenario to output DALY with
        inputs_aggregated (bool, optional): Flag to indicate that inputs are aggregated. Defaults to True.
        save_number (bool, optional): Flag to save numbers. Defaults to True.
        provenance (bool, optional): Flag to indicate saving with Provenance. Defaults to False.
        holds (List, optional): Optional Holds to pass to submit_cluster_job_request. Defaults to [1].
        extra_flags (List, optional): Custom arguments to pass to submit_cluster_job_request. Defaults to [].
        memory (int, optional): Memory to pass to submit_cluster_job_request. Defaults to 100.
        queue (str, optional): Queue to pass to submit_cluster_job_request. Defaults to "all".
        runtime (str, optional): Runtime to pass to submit_cluster_job_request. Defaults to "04:00:00".
        job_name (str, optional): Job name to pass to submit_cluster_job_request. Defaults to "calculate_dalys".

    Returns:
        int: Job ID for submitted request.
    """
    dalys_args = {
        "gbd_round_id": 6,
        "years": years,
        "draws": draws,
        "versions": versions,
        "output_scenario": output_scenario,
        "inputs_aggregated": inputs_aggregated,
        "save_number": save_number,
        "provenance": provenance,
    }

    command_args = args_from_ctx(dalys_args)

    executable_path = str(Path(sys.executable).parent / "fhs_pipeline_dalys_console")
    summarize_call = f"{executable_path} parallelize-by-cause {command_args}"

    return submit_cluster_job_request(
        job_name=job_name,
        memory=memory,
        cores=1,
        command_to_run=summarize_call,
        runtime=runtime,
        cluster_queue=f"{queue}.q",
        job_holds=holds,
        extra_flags=extra_flags,
    )
