import pandas as pd

from get_draws.api import get_draws

from cod_prep.claude.configurator import Configurator

CONF = Configurator("standard")
DEM_COLS = ["cause_id", "location_id", "sex_id", "year_id", "age_group_id"]
DRAW_COLS = ["draw_" + str(x) for x in range(0, CONF.get_id("draw_num"))]


def calculate_relative_risk(cfr_exposed, cfr_unexposed):
    """Relative risk of case fatality ratios (CFR)."""
    df = cfr_exposed.merge(cfr_unexposed, how="left", on=DEM_COLS)
    df[DRAW_COLS] = df.filter(regex="draw_.*_x").rename(columns=lambda c: c[:-2]) / df.filter(
        regex="draw_.*_y"
    ).rename(columns=lambda c: c[:-2])
    return df[DEM_COLS + DRAW_COLS]


def calculate_paf(rr, incidence):
    """Population attributable fraction.

    The amount of disease incidence that is attributed to a specific risk
    factor or the amount of potential decrease in disease incidence that can
    be expected if that risk factor were completely eliminated
    https://www.who.int/quantifying_ehimpacts/faqs/en/

    PAF = p * (rr-1) / (1 + p * (rr-1))
    p: proportion of the population exposed
    rr: relative risk
    """
    incidence = add_population(
        incidence,
        force_rerun=False,
        block_rerun=True,
        population_run_id=CONF.get_id("pop_run"),
    )
    # are there draws of population?
    p = incidence.filter(regex="draw_") / incidence["mean"]
    # p[DRAW_COLS] = p[DRAW_COLS] / 100000   why is that here???
    df = p.merge(rr, how="left", on=DEM_COLS)
    for draw in DRAW_COLS:
        df[draw] = (df[f"{draw}_x"] * (df[f"{draw}_y"] - 1)) / (
            1 + (df[f"{draw}_x"] * (df[f"{draw}_y"] - 1))
        )
    return df[DEM_COLS + DRAW_COLS]


def calculate_attr_burden(paf):
    """Attributable burden calculation for deaths & YLLs."""
    cc_df = get_draws(
        gbd_id_type="cause_id",
        gbd_id=cause_id,
        source="codcorrect",
        metric_id=1,
        measure_id=[1, 4],
        location_id=list(paf.location_id.unique()),
        year_id=year_id,
        gbd_round_id=CONF.get_id("gbd_round"),
        num_workers=2,
        version_id=CONF.get_id("codcorrect_version"),
        decomp_step=CONF.get_id("decomp_step"),
    )
    df = cc_df.merge(paf, how="inner", on=DEM_COLS)
    df[DRAW_COLS] = df.filter(regex="draw_.*_x").rename(columns=lambda c: c[:-2]) * df.filter(
        regex="draw_.*_y"
    ).rename(columns=lambda c: c[:-2])
    return df[DEM_COLS + DRAW_COLS + ["measure_id"]]


def main():
    cfr_exposed = pd.read_csv(f"{base_dir}/{model_run}/{year_id}/{cause_id}_prop.csv")
    cfr_unexposed = pd.read_csv(f"{base_dir}/{model_run}/{year_id}/{cause_id}_prop.csv")
    rr = calculate_relative_risk(cfr_exposed, cfr_unexposed)
    incidence = pd.read_csv(f"{base_dir}/{model_run}/{year_id}/{cause_id}_aggregate.csv")
    paf = calculate_paf(rr, incidence)
    attr_burden = calculate_attr_burden(paf)
    attr_burden.to_csv()
