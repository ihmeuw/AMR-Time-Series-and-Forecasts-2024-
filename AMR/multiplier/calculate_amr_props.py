import pandas as pd
import numpy as np
import sys
import os
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from cod_prep.utils import (
    print_log_message, get_norway_subnational_mapping,
    report_if_merge_fail)
from pathlib import Path
from cod_prep.downloaders import (
    get_current_location_hierarchy, add_location_metadata,
    add_population, get_current_cause_hierarchy)
import warnings
from get_draws.api import get_draws
from db_queries import get_outputs
from db_queries import get_model_results
from amr_prep.utils.misc import get_prepped_csbg_universe
from cod_prep.utils.misc import report_duplicates

CONF = Configurator()


class AmrPropCalculator(object):

    def __init__(self, burden_type, year_id, pathogen, aggregate_fraction_of_resistance=False):
        self.burden_type = burden_type
        self.year_id = year_id
        self.pathogen = pathogen
        self.aggregate_fraction_of_resistance = aggregate_fraction_of_resistance
        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.dem_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']
        self.draw_cols = ["draw_" + str(x) for x in range(0, self.conf.get_id('draw_num'))]
        self.extra_draw_cols = ["draw_" + str(x) for x in range(100, 1000)]

    def subset_or_exceptions(self, df, step, abx_class):
        df = df.copy()
        if self.pathogen in [
            'enteropathogenic_escherichia_coli',
            'enterotoxigenic_escherichia_coli'
        ]:
            return df.query(
                f"pathogen == 'escherichia_coli' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen)
        elif (
            (self.pathogen == 'salmonella_paratyphi') & (abx_class != 'mdr') & (step == 'rr')
        ):
            return df.query(
                f"pathogen == 'salmonella_typhi' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen)
        elif (
            (self.pathogen in ['salmonella_paratyphi', 'salmonella_typhi']) & (
                abx_class == 'mdr') & (step == 'rr')
        ):
            return df.query(
                f"pathogen == 'salmonella_typhi' & abx_class == 'fluoroquinolone'"
            ).assign(pathogen=self.pathogen, abx_class=abx_class)
        elif (
            (self.pathogen == 'shigella_spp') & (abx_class == 'fluoroquinolone') & (
                step == 'rr') & (self.burden_type == 'fatal')
        ):
            return df.query(
                f"pathogen == 'salmonella_typhi' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen)
        elif (
            (self.pathogen == 'salmonella_ints_non_typhi') &
            (step == 'rr')
        ):
            return df.query(
                f"pathogen == 'salmonella_typhi' & abx_class == 'fluoroquinolone'"
            ).assign(pathogen=self.pathogen, abx_class=abx_class)
        elif (
            (self.pathogen == 'streptococcus_group_a') & (step == 'rr')
        ):
            return df.query(
                f"pathogen == 'group_a_strep' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen, abx_class=abx_class)
        elif (
            (self.pathogen == 'streptococcus_group_b') & (step == 'rr')
        ):
            return df.query(
                f"pathogen == 'group_b_strep' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen, abx_class=abx_class)
        else:
            return df.query(
                f"pathogen == '{self.pathogen}' & abx_class == '{abx_class}'"
            )

    def get_frac_resist_ids(self, abx_class):
        if not (
            (self.pathogen == 'mycobacterium_tuberculosis') & (
                abx_class in ['mdr', 'xdr'])
        ):
            bug_drug_universe = pd.read_csv("FILEPATH", encoding='latin1')
            bug_drug = self.subset_or_exceptions(bug_drug_universe, 'pr', abx_class)
            assert len(bug_drug) == 1,\
                f"{self.pathogen}/{abx_class} is not a valid bug/drug"
            return [int(bug_drug['resistance_run_id'].unique()[0]), int(bug_drug['modelable_entity_id'].unique()[0])]
        else:
            return None, None

    def get_stgpr_run_params(self, run_id):
        stgpr_run_params = pd.read_csv("FILEPATH")
        assert len(stgpr_run_params) == 1
        return stgpr_run_params.loc[0].to_dict()

    def get_stgpr_fraction_of_resistance(self, abx_class, run_id, modeleable_entity_id):
        if self.pathogen == 'mycobacterium_tuberculosis':
            tb_path = "FILEPATH"
            df = pd.read_csv(tb_path + '.csv')
            df = df.loc[df['abx_class'] == abx_class, ]
        else:
            run_params = self.get_stgpr_run_params(run_id)
            st_gpr_draws = Path("FILEPATH")
            lh = get_current_location_hierarchy(
                location_set_version_id=self.conf.get_id("location_set_version_2023"),
                **self.cache_kwargs
            ).query("is_estimate == 1").assign(bad_loc=0)

            if run_params['gbd_round_id'] == 7:
                nor = get_norway_subnational_mapping().set_index(
                    'location_id_old')['location_id_new'].to_dict()
                lh['location_id'] = lh['location_id'].replace(nor)
            elif run_params['release_id'] == 16:
                lh = lh.loc[~(lh['location_id'].isin([60908, 94364])), ]
            dfs = []
            for location_id in lh.location_id.unique().tolist():
                loc_file = "FILEPATH"
                if not loc_file.exists():
                    lh.loc[lh.location_id == location_id, 'bad_loc'] = 1
                else:
                    loc_df = pd.read_csv(loc_file)
                    loc_df = loc_df.query(f"year_id == {self.year_id}")
                    dfs.append(loc_df)
            df = pd.concat(dfs)

            stgpr_pe = get_model_results(
                gbd_team='stgpr',
                gbd_id=modeleable_entity_id,
                model_version_id=run_id,
                release_id=run_params['release_id'],
                year_id=self.year_id,
            )[self.dem_cols + ['mean']]
            stgpr_pe.rename(columns={'mean': 'point_estimate'}, inplace=True)
            df = df.merge(stgpr_pe, on=self.dem_cols, how='left', validate='one_to_one')    

            df.loc[df['location_id'] == 95069, 'location_id'] = 44858

            if lh.bad_loc.sum() > 0:
                raise AssertionError(
                    f"The following locations are missing:\n "
                    f"{lh.loc[lh.bad_loc == 1, 'location_name'].unique().tolist()}"
                )
            if run_params['gbd_round_id'] == 7:
                nor = get_norway_subnational_mapping().rename(
                    columns={'location_id_new': 'location_id'}
                )
                df = df.merge(nor, how='left', on=['location_id'])
                df.loc[df.location_id_old.notnull(), 'location_id'] = df['location_id_old']
            df = df[self.dem_cols + self.draw_cols + ['point_estimate']]

            assert run_params['draws'] == self.conf.get_id('draw_num')
            assert run_params['data_transform'] == 'logit'
            df[self.draw_cols + ['point_estimate']] = np.clip(df[self.draw_cols + ['point_estimate']], 0, 1, None)
            assert (df[self.draw_cols + ['point_estimate']] >= 0).values.all()
            assert (df[self.draw_cols + ['point_estimate']] <= 1).values.all()

        if self.aggregate_fraction_of_resistance:
            print_log_message("Aggregating to super region")
            df = add_location_metadata(
                df, 'super_region_id',
                location_set_version_id=self.conf.get_id('location_set_version'),
                **self.cache_kwargs
            )
            report_if_merge_fail(df, 'super_region_id', 'location_id')
            df = add_population(
                df, pop_run_id=self.conf.get_id('pop_run'), **self.cache_kwargs
            )
            report_if_merge_fail(df, 'population', 'location_id')
            group_cols = ['super_region_id', 'age_group_id', 'sex_id', 'year_id']
            df['weight'] = df['population'] / df.groupby(group_cols).population.transform(sum)
            df[self.draw_cols + ['point_estimate']] = pd.DataFrame(
                df[self.draw_cols + ['point_estimate']].to_numpy() * df[['weight']].to_numpy(),
                index=df.index
            )
            df[self.draw_cols] = df.groupby(group_cols)[self.draw_cols + ['point_estimate']].transform(sum)
            assert (df[self.draw_cols + ['point_estimate']] <= 1).values.all()
            df = df.drop(['super_region_id', 'population', 'weight'], axis='columns')
        df = df.assign(pathogen=self.pathogen, abx_class=abx_class)
        return df

    def get_resistance_profiles(self):
        lh = get_current_location_hierarchy(
            location_set_version_id=self.conf.get_id("location_set_version"),
            **self.cache_kwargs
        ).query("level == 3").assign(bad_loc=0)
        profile_path = self.conf.get_resource("resistance_profile")
        loc_dfs = []
        for location_id in lh.location_id.unique().tolist():
            if self.pathogen in [
                'enteropathogenic_escherichia_coli',
                'enterotoxigenic_escherichia_coli'
            ]:
                pathogen_profile = 'escherichia_coli'
            elif self.pathogen == 'streptococcus_group_b':
                pathogen_profile = 'group_b_strep'
            else:
                pathogen_profile = self.pathogen
            loc_file = Path(profile_path.format(
                pathogen=pathogen_profile, location_id=location_id
            ))
            if not loc_file.exists():
                lh.loc[lh.location_id == location_id, 'bad_loc'] = 1
            else:
                loc_dfs.append(pd.read_csv(loc_file).query(f"year_id == {self.year_id}"))
        if lh.bad_loc.sum() > 0:
            raise AssertionError(
                f"The following locations are missing:\n "
                f"{lh.loc[lh.bad_loc == 1, 'location_name'].unique().tolist()}"
            )
        df = pd.concat(loc_dfs, sort=False)

        df.loc[df['pathogen'] == 'group_b_strep', 'pathogen'] = 'streptococcus_group_b'
        assert {'age_group_id', 'sex_id'}.isdisjoint(df)
        df = df.assign(age_group_id=22, sex_id=3)
        df['point_estimate'] = df[self.draw_cols].mean(1)
        return df

    def get_relative_risk_of_death(self, abx_classes):
        df = pd.read_csv("FILEPATH")

        df = pd.concat([
            self.subset_or_exceptions(df, 'rr', abx_class)
            for abx_class in abx_classes
        ], sort=False)
        df = df.rename(
            lambda x: 'draw_' + str(int(x.replace("X", "")) - 1) if 'X' in x else x,
            axis='columns'
        )
        df = df[['infectious_syndrome', 'pathogen', 'abx_class'] + self.draw_cols + ['point_estimate']]
        df[self.draw_cols + ['point_estimate']] = np.exp(df[self.draw_cols + ['point_estimate']])

        if self.pathogen == 'mycobacterium_tuberculosis':
            base_dir = "FILEPATH"
            df = pd.read_csv("FILEPATH")
            df['year_id'] = self.year_id
            df['infectious_syndrome'] = 'L2_tb'

        assert len(df['abx_class'].unique()) == len(abx_classes)
        assert (df[self.draw_cols + ['point_estimate']] >= 0).values.all()
        if (df[self.draw_cols + ['point_estimate']] < 1).values.any():
            warnings.warn("Some draws show a protective relative risk")
        return df

    def get_relative_risk_of_los(self, abx_classes):
        all_los = pd.read_csv("FILEPATH")
        dfs = []
        for abx_class in abx_classes:
            if (self.pathogen == 'mycobacterium_tuberculosis') & (
                abx_class in ['mdr', 'xdr']
            ):
                lh = get_current_location_hierarchy(**self.cache_kwargs)
                locs = lh.loc[lh['level'] == 3, 'location_id'].unique().tolist()
                if abx_class == 'mdr':
                    exposed = 946
                    unexposed = 934
                elif abx_class == 'xdr':
                    exposed = 947
                    unexposed = 934
                draw_kwargs = {
                    'gbd_id_type': 'cause_id',
                    'gbd_id': [exposed, unexposed],
                    'source': 'como',
                    'version_id': self.conf.get_id("como_version"),
                    'downsample': True,
                    'n_draws': CONF.get_id('draw_num'),
                }
                output_kwargs = {
                    'topic': 'cause',
                    'compare_version_id': CONF.get_id('compare_version'),
                    'cause_id': [exposed, unexposed],
                }
                kwargs = {
                    'measure_id': [3, 6],
                    'location_id': 1,
                    'year_id': self.year_id,
                    'age_group_id': 22,
                    'sex_id': 3,
                    'metric_id': 3,
                    'release_id': self.conf.get_id("release"),
                }
                print_log_message("Pulling COMO results for TB LOS RR")
                df = get_draws(**{**draw_kwargs, **kwargs})
                pe = get_outputs(**{**output_kwargs, **kwargs})
                pe.rename(columns={'val': 'point_estimate'}, inplace=True)
                pe_keep_cols = ['location_id', 'age_group_id', 'year_id', 'measure_id', 'metric_id', 'cause_id']
                pe = pe[pe_keep_cols + ['point_estimate']]
                df = df.merge(pe, on=pe_keep_cols, how='left', validate='one_to_one')
                yld_exposed = df.query(
                    f"measure_id == 3 & cause_id == {exposed}").reset_index()[self.draw_cols + ['point_estimate']]
                cases_exposed = df.query(
                    f"measure_id == 6 & cause_id == {exposed}").reset_index()[self.draw_cols + ['point_estimate']]
                yld_unexposed = df.query(
                    f"measure_id == 3 & cause_id == {unexposed}").reset_index()[self.draw_cols + ['point_estimate']]
                cases_unexposed = df.query(
                    f"measure_id == 6 & cause_id == {unexposed}").reset_index()[self.draw_cols + ['point_estimate']]
                df = yld_exposed.div(cases_exposed).div(yld_unexposed.div(cases_unexposed))
                df = df.assign(pathogen=self.pathogen, abx_class=abx_class)
                if self.year_id == 1990:
                    df = df.fillna(0)
            else:
                df = self.subset_or_exceptions(all_los, 'rr', abx_class)
                if len(df) == 0:
                    df = all_los_prior.query(f"abx_class == '{abx_class}'")
                    df = df.assign(pathogen=self.pathogen)
                assert len(df) == 1
                df = df.rename(
                    lambda x: 'draw_' + str(int(x) - 1) if x not in ['point_estimate', 'pathogen', 'abx_class'] else x,
                    axis='columns'
                )
                df[self.draw_cols + ['point_estimate']] = np.exp(df[self.draw_cols + ['point_estimate']])
            dfs.append(df)

        df = pd.concat(dfs, sort=False)
        assert len(df) == len(abx_classes)
        assert (df[self.draw_cols + ['point_estimate']] >= 0).values.all()
        if (df[self.draw_cols + ['point_estimate']] < 1).values.any():
            warnings.warn("Some draws show a protective relative risk")

        if 'infectious_syndrome' not in df.columns:
            df['infectious_syndrome'] = 'all'

        return df

    def process_inputs(self, resist_cases, rr):
        rr = pd.merge(
            resist_cases[
                ['pathogen', 'abx_set', 'combinatoric_id', 'abx_class']
            ].drop_duplicates(), rr, how='outer',
            on=['pathogen', 'abx_class'],
            indicator=True
        )
        assert (rr._merge == 'both').all()
        rr = rr.drop('_merge', axis='columns')
        rr_group_cols = ['infectious_syndrome', 'pathogen', 'abx_set', 'combinatoric_id']
        rr_prof = rr.groupby(rr_group_cols, as_index=False)[self.draw_cols + ['point_estimate']].max()

        resist_cases = resist_cases.drop('abx_class', axis=1)\
            .drop_duplicates().reset_index(drop=True)
        resist_unique_cols = [
            'pathogen', 'abx_set', 'combinatoric_id',
            'location_id', 'year_id', 'age_group_id', 'sex_id'
        ]
        rr_unique_cols = ['infectious_syndrome', 'pathogen', 'abx_set', 'combinatoric_id']
        if (pathogen == 'mycobacterium_tuberculosis'):
            resist_unique_cols += ['measure_id']

        report_duplicates(resist_cases, resist_unique_cols)
        report_duplicates(rr_prof, rr_unique_cols)
        resist_cases = resist_cases.set_index(resist_unique_cols)
        rr_prof = rr_prof.set_index(rr_unique_cols)

        if pathogen != 'mycobacterium_tuberculosis':
            reall_props = rr.copy().reset_index(drop=True)

            beta_lactams = [
                'penicillin', 'aminopenicillin', 'beta_lactamase_inhibitor',
                'third_gen_ceph', 'fourth_gen_ceph',
                'anti_pseudomonal_penicillin', 'carbapenem'
            ]
            beta_lactams_dict = dict(zip(beta_lactams, list(range(1, 8))))
            reall_props['beta_lactams_rank'] = reall_props['abx_class'].map(beta_lactams_dict)
            reall_props['combo_max_bl'] = reall_props.groupby(
                ['combinatoric_id'], as_index=False).beta_lactams_rank.transform('max')
            reall_props = reall_props.loc[
                (reall_props['beta_lactams_rank'].isna()) | (
                    reall_props['beta_lactams_rank'] == reall_props['combo_max_bl']),
            ]
            reall_props.drop(columns=['beta_lactams_rank', 'combo_max_bl'], inplace=True)
            reall_props[self.draw_cols + ['point_estimate']] = reall_props[self.draw_cols + ['point_estimate']] - 1
            all_protective = reall_props.groupby(
                ['infectious_syndrome', 'pathogen', 'abx_set', 'combinatoric_id']
            )[self.draw_cols + ['point_estimate']].transform(lambda x: (x < 0).all())
            mask = (~all_protective & (reall_props[self.draw_cols + ['point_estimate']] < 0)).to_numpy()
            props = reall_props[self.draw_cols + ['point_estimate']].to_numpy()
            props[mask] = 0
            mask2 = all_protective.to_numpy()
            props[mask2] = 0.001
            reall_props[self.draw_cols + ['point_estimate']] = pd.DataFrame(props, index=reall_props.index)
            reall_props[self.draw_cols + ['point_estimate']] = reall_props[self.draw_cols + ['point_estimate']] / reall_props.groupby(
                ['infectious_syndrome', 'pathogen', 'abx_set', 'combinatoric_id']
            )[self.draw_cols + ['point_estimate']].transform(sum)
            assert (reall_props[self.draw_cols + ['point_estimate']] >= 0).values.all()
            assert (reall_props[self.draw_cols + ['point_estimate']] <= 1).values.all()
            assert (reall_props != np.inf).values.all()
            reall_props = reall_props.set_index(
                ['pathogen', 'abx_set', 'combinatoric_id', 'abx_class', 'infectious_syndrome']
            )
            self.reall_props = reall_props
        self.resist_unique_cols = resist_unique_cols
        self.resist_cases = resist_cases
        self.rr_prof = rr_prof
        

    def calculate_props(self, counterfactual, measure_id):
        print_log_message(
            f"Calculating props for {counterfactual}, measure {measure_id}"
        )
        assert counterfactual in ['no_infection', 'susceptible_infection']
        assert measure_id in [1, 3, 4, 6]
 
        if measure_id in [1, 3, 4]:
            def reord(df, addition_col=['infectious_syndrome']):
                return df.reorder_levels(
                    [c for c in self.resist_unique_cols if c in df.index.names] + addition_col
                )

            group_cols = [
                'pathogen', 'abx_set', 'location_id', 'year_id',
                'age_group_id', 'sex_id'
            ]
            if counterfactual == 'no_infection':
                if (self.pathogen == 'mycobacterium_tuberculosis') & (self.burden_type == 'fatal'):
                    df = self.resist_cases.reset_index()
                    df = df.loc[df['measure_id'] == 1, ].set_index(self.resist_unique_cols)
                else:
                    if (self.pathogen == 'mycobacterium_tuberculosis') & (self.burden_type == 'nonfatal'):
                        tb_rc = self.resist_cases.reset_index()
                        self.resist_cases = tb_rc.loc[tb_rc['measure_id'] == 6, ].set_index(self.resist_unique_cols)
                    total_resist = reord(
                        self.resist_cases.groupby(level=group_cols)[self.draw_cols + ['point_estimate']].sum(),
                        []
                    )
                    total_relative = reord(self.resist_cases.mul(self.rr_prof).groupby(
                        level=group_cols + ['infectious_syndrome'])[self.draw_cols + ['point_estimate']].sum(),
                    )
                    df = reord(
                        reord(self.resist_cases.mul(self.rr_prof)).div(
                            reord(total_relative.add(1 - total_resist))
                        )
                    )
            elif counterfactual == 'susceptible_infection':
                if pathogen == 'mycobacterium_tuberculosis':
                    incd_res = self.resist_cases.reset_index()
                    incd_res = incd_res.loc[incd_res['measure_id'] == 6, ].set_index(self.resist_unique_cols)
                    excess_risk = reord(incd_res.mul(self.rr_prof - 1))
                else:
                    excess_risk = reord(self.resist_cases.mul(self.rr_prof - 1))
                total_excess_risk = reord(
                    excess_risk.groupby(level=group_cols + ['infectious_syndrome'])[self.draw_cols + ['point_estimate']].sum())
                df = reord(excess_risk.div(1 + total_excess_risk))
        elif measure_id == 6:
            if counterfactual == 'no_infection':
                df = self.resist_cases.copy()
                if self.pathogen == 'mycobacterium_tuberculosis':
                    df = df.reset_index()
                    df = df.loc[df['measure_id'] == 6, ].set_index(self.resist_unique_cols)
            elif counterfactual == 'susceptible_infection':
                df = self.resist_cases.copy()
                if self.pathogen == 'mycobacterium_tuberculosis':
                    df = df.reset_index()
                    df = df.loc[df['measure_id'] == 6, ].set_index(self.resist_unique_cols)
                df[self.draw_cols + ['point_estimate']] = 0

        assert (df[self.draw_cols + ['point_estimate']] <= 1).values.all()
        if counterfactual == 'no_infection':
            assert (df[self.draw_cols + ['point_estimate']] >= 0).values.all()
        elif counterfactual == 'susceptible_infection':
            assert (df[self.draw_cols] >= -1).values.all()
        assert df.notnull().values.all()
        assert (df != np.inf).values.all()

        if pathogen != 'mycobacterium_tuberculosis':
            print_log_message("Reallocating props back to drugs")
            df = df.mul(self.reall_props)
            df = df.groupby(
                level=['infectious_syndrome', 'pathogen', 'abx_set', 'abx_class'] + self.dem_cols
            )[self.draw_cols + ['point_estimate']].sum()
            df = df.reset_index(drop=False)
        else:
            df = df.reset_index(drop=False)
            if 'infectious_syndrome' not in df.columns:
                df['infectious_syndrome'] = 'L2_tb'
            df.loc[df['abx_set'] == 'mdr_only', 'abx_class'] = 'mdr'
            df.loc[df['abx_set'] == 'xdr_only', 'abx_class'] = 'xdr'

        assert df.notnull().values.all()
        df = df[self.dem_cols + ['infectious_syndrome', 'pathogen', 'abx_set', 'abx_class'] + self.draw_cols + ['point_estimate']]
        df = df.assign(counterfactual=counterfactual, measure_id=measure_id)
        return df

    def reconcile_counterfactuals(self, df):

        susc = df.loc[df.counterfactual == 'susceptible_infection']
        df = df.loc[df.counterfactual != 'susceptible_infection']
        assert (susc.abx_set == 'all').all()
        assert not susc.abx_class.isin(['all_susceptible', 'all_resistant']).any()

        merge_cols = self.dem_cols + ['pathogen', 'abx_class', 'measure_id', 'infectious_syndrome']
        no_inf = pd.merge(
            df[merge_cols + self.draw_cols + ['point_estimate']],
            susc[merge_cols], how='right', on=merge_cols,
            validate='one_to_one'
        )
        susc = susc.drop(['abx_set', 'counterfactual'], axis='columns')
        susc = susc.set_index(merge_cols)
        no_inf = no_inf.set_index(merge_cols)
        susc = susc.mask(no_inf.sub(susc) < 0, other=no_inf)
        susc = susc.reset_index().assign(
            abx_set='all', counterfactual='susceptible_infection'
        )[df.columns.tolist()]
        df = df.append(susc)
        assert df.notnull().values.all()
        return df

    def run_calculator(self):
        print_log_message("Loading inputs")
        print_log_message("Getting fraction of resistance in cases")
        csbg = get_prepped_csbg_universe(self.burden_type, add_parents=False)
        if self.pathogen == 'acinetobacter_baumanii':
            name_fix = 'acinetobacter_baumannii'
        else:
            name_fix = self.pathogen
        abxs = csbg.loc[
            csbg.abx_class.notnull() & (csbg.pathogen == name_fix), 'abx_class'
        ].unique().tolist()

        resist_cases = pd.concat([
            self.get_stgpr_fraction_of_resistance(
                abx_class, self.get_frac_resist_ids(abx_class)[0], self.get_frac_resist_ids(abx_class)[1]
            ) for abx_class in abxs
        ], sort=False).assign(
            abx_set=lambda d: d['abx_class'] + '_only',
            combinatoric_id=lambda d: d['pathogen'] + '-' + d['abx_class'],
        )
        if self.pathogen != 'mycobacterium_tuberculosis':
            if len(abxs) > 1:
                print_log_message("Getting fractions for resistance profiles")
                resist_prof = self.get_resistance_profiles()
                assert set(resist_prof.abx_class) == set(abxs)
                assert len(set(resist_prof.combinatoric_id)) == 2 ** len(abxs)
                assert np.allclose(
                    resist_prof.groupby(
                        ['pathogen', 'location_id', 'age_group_id', 'sex_id',
                         'year_id', 'abx_class']
                    )[self.draw_cols + ['point_estimate']].sum(), 1)
                resist_prof = resist_prof.query("resistant == 1").drop(
                    'resistant', axis='columns'
                )
                if self.pathogen in [
                    'enteropathogenic_escherichia_coli',
                    'enterotoxigenic_escherichia_coli'
                ]:
                    resist_prof['pathogen'] = self.pathogen
                    resist_prof['combinatoric_id'] = resist_prof['combinatoric_id'].str.replace(
                        'escherichia_coli', self.pathogen
                    )
                resist_cases = resist_cases.append(resist_prof.assign(abx_set='all'))
            else:
                resist_cases = resist_cases.append(resist_cases.assign(abx_set='all'))

        if self.burden_type == 'fatal':
            print_log_message("Getting relative risk of death")
            rr = self.get_relative_risk_of_death(abxs)
            measures = [1, 4]
        elif self.burden_type == 'nonfatal':
            print_log_message("Getting relative risk of LOS")
            rr = self.get_relative_risk_of_los(abxs)
            measures = [3, 6]

        print_log_message("Prepping inputs")
        self.process_inputs(resist_cases, rr)

        print_log_message("Calculating AMR props for 2 counterfactuals")
        df = pd.concat([
            self.calculate_props(counterfactual, measure_id)
            for counterfactual in ['no_infection', 'susceptible_infection']
            for measure_id in measures
        ], sort=False)
        non_draw_cols = list(set(df.columns.values) - set(self.draw_cols + ['point_estimate']))
        non_draw_cols.remove('abx_class')
        
        df = df.loc[~((df['measure_id'] == 6) & (df['counterfactual'] == 'susceptible_infection')), ]

        if (self.pathogen == 'mycobacterium_tuberculosis'):
            df['abx_set'] = 'all'
            no_inf_tb = df.loc[df['counterfactual'] == 'no_infection', ]
            AR_noinf_tb = no_inf_tb.groupby(non_draw_cols)[self.draw_cols + ['point_estimate']].sum()
            susc_tb = (1 - AR_noinf_tb)
            susc_tb.reset_index(inplace=True)
            AR_noinf_tb.reset_index(inplace=True)
            susc_tb['abx_class'] = 'all_susceptible'
            AR_noinf_tb['abx_class'] = 'all_resistant'
            df = pd.concat([df, AR_noinf_tb, susc_tb])
        else:
            no_ddp_with_resis = df.loc[
                (df['counterfactual'] == 'no_infection') & (
                    df['abx_set'] != 'all'), ]
            ddp_with_resis = df.loc[
                (df['counterfactual'] == 'no_infection') & (
                    df['abx_set'] == 'all'), ]
            df = df.loc[(df['counterfactual'] != 'no_infection') & (df['abx_set'] == 'all'), ]
            ddp_with_resis = ddp_with_resis.groupby(
                non_draw_cols, as_index=False)[self.draw_cols + ['point_estimate']].sum()
            assert (ddp_with_resis[self.draw_cols + ['point_estimate']] <= 1).values.all()
            ddp_with_resis.set_index(non_draw_cols, inplace=True)
            ddp_susc = (1 - ddp_with_resis)
            ddp_susc.reset_index(inplace=True)
            ddp_with_resis.reset_index(inplace=True)
            ddp_susc['abx_class'] = 'all_susceptible'
            ddp_with_resis['abx_class'] = 'all_resistant'
            df = pd.concat([df, no_ddp_with_resis, ddp_with_resis, ddp_susc])

        df = self.reconcile_counterfactuals(df)
        print_log_message("Saving results")
        AmrResult(
            process='calculate_amr_props',
            burden_type=self.burden_type,
            year_id=self.year_id,
            pathogen=self.pathogen,
        ).write_results(df)


if __name__ == '__main__':
    burden_type = str(sys.argv[1])
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        base_dir = "FILEPATH"
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        pathogen = str(task_row['pathogen'])
    else:
        year_id = int(sys.argv[2])
        pathogen = str(sys.argv[3])
    calc = AmrPropCalculator(
        burden_type=burden_type,
        year_id=year_id,
        pathogen=pathogen
    )
    calc.run_calculator()
