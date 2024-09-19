#imports
import pandas as pd
import numpy as np
from pathlib import Path
from db_tools import ezfuncs
from cod_prep.downloaders import *
from cod_prep.utils import *
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import get_amr_results
from amr_prep.utils.amr_io import AmrResult
from db_queries import get_population 
from mcod_prep.utils.mcause_io import get_mcause_results
import itertools
import string

CONF = Configurator()
LSV_ID = CONF.get_id('location_set_version')
CSV_ID = CONF.get_id('cause_set_version')

lh = get_current_location_hierarchy(location_set_version_id=LSV_ID)
ch = get_current_cause_hierarchy(cause_set_version_id=CSV_ID)

#df =pd.read_parquet({FILEPATH}) 
all_age_underlying_comp = df.loc[df["age_group_id"]==22]

ages = get_cod_ages()
cod_ages = ages.age_group_id.unique().tolist()
df=df.loc[df["age_group_id"].isin(cod_ages)]

class das_gupta_decomp:
    def __init__(self, N:int):
        
        # Checking if N exceeds 26
        if N > 26:
            raise ValueError(f"Cannon create decomposition equation for {N} factors. Up to 26 factors are currently supported.")

        # Setting up self-variables
        self.__N = N
        self.__R1_factors = list(string.ascii_uppercase)[:N]
        self.__R2_factors = list(string.ascii_lowercase)[:N]
        self._effect_eqs = self.get_equations()
    
    def __str__(self):
        return f"<Das Gupta Decomposition Object of {self._N} Factors>"
    
    
    def get_equations(self):
        # From number of factors, determining number of terms in Das Gupta decomposition equation
        terms = self.__N // 2 + self.__N % 2

        # Initializing dictionary containing equations
        effect_eqs = {}

        # For each decomp factor...
        for idx, factor in enumerate(self.__R1_factors):

            # Getting factors for each rate, not including the current factor being analyzed
            R1_facts = [elem for elem in self.__R1_factors if elem!=factor]
            R2_facts = [elem for elem in self.__R2_factors if elem!=factor.lower()]

            # Getting valid combinations of factors for both rates
            # Valid combinations are any combinations of R1 and R2 elements
            # Of length N-1 where all elements in that combination are non-duplicate letters
            # ex. 'B' and 'b' are considered duplicate letters
            combos = [combo for combo in itertools.combinations(R1_facts+R2_facts, self.__N-1) if len(set([x.lower() for x in combo]))==self.__N-1]

            # Adding a count of the number of capital factors (need to change name)
            n_caps = {f"({'*'.join(combo)})":sum(1 for c in combo if c.isupper()) for combo in combos}

            # Creating list to contain equation terms
            term_list = []

            # For each term in equation...
            for t in range(terms):

                # Getting term numerator
                numerator = '+'.join([k for k,v in n_caps.items() if v==(self.__N-(self.__N-t)) or v==(self.__N-(t+1))])

                # Getting term denominator
                denominator = str(int(np.prod([x for x in range(self.__N-t,self.__N+1)]) / max(1,t)))

                # Adding term to term list
                term_list += [f"(({numerator})/({denominator}))"]

            # Creating equation and adding to equation dictionary
            equation = f"({'+'.join(term_list)})*({factor.lower()}-{factor})"
            #equation = f"({'+'.join(term_list)})"
            effect_eqs[f"{factor}_effect"] = equation

        return effect_eqs

    def decompose(self, base_df: pd.DataFrame, comp_df: pd.DataFrame, col_names:list=[], base_col_names:list=[], comp_col_names:list=[]):
        
        # Raising error if rows are not equal
        if len(base_df) != len(comp_df):
            raise RuntimeError("Dataframes must have the same number of rows. Decompositions are calculated for each row pair.")
        
        # Raising error if too many parameters are provided
        if not (col_names or base_col_names or comp_col_names):
            raise ValueError("Must provide col_names OR base_col_names and comp_col_names.")
        elif col_names and (base_col_names or comp_col_names):
            raise ValueError("Must only provide either col_names or base_col_names and comp_col_names.")
        elif not col_names and not (base_col_names and comp_col_names):
            raise ValueError("Must only provide both base_col_names and comp_col_names if one is provided.")
        
        # Creating name_dict object to map from base_df to comp_df columns
        if base_col_names and comp_col_names:
            name_dict = dict(zip(base_col_names,comp_col_names))
        else:
            name_dict = dict(zip(col_names,col_names))
        
        # Raising error if the wrong number of columns were provided
        # Also throws error if duplicate column names are provided
        if len(name_dict) != self.__N:
            raise ValueError(f"Expecting to decompose {self.__N} columns but {len(name_dict)} were provided")
        
        # Raising error if any columns provided are not present in either dataframe
        if not all([c in base_df.columns for c in name_dict.keys()]):
            raise KeyError("Some column names provided are not present in base_df.")
        elif not all([c  in comp_df.columns for c in name_dict.values()]):
            raise KeyError("Some column names provided are not present in comp_df.")
            
        # Mapping variable names to corresponding columns
        var_to_name = {var:col for var,col in zip(self.__R1_factors+self.__R2_factors, list(name_dict.keys())+list(name_dict.values()))}
        var_to_col = {var:col for var,col in zip(self.__R1_factors+self.__R2_factors, [f"base_df['{i}'].values" for i in name_dict.keys()]+[f"comp_df['{i}'].values" for i in name_dict.values()])}
        
        
        # Getting translated equations
        trans_effect_eqs = {f"{self._translate_string(eff[0], var_to_name)}_effect":self._translate_string(eq, var_to_col) for eff, eq in self._effect_eqs.items()}
        
        decomp_df = pd.DataFrame()
        for f, eq in trans_effect_eqs.items():
            decomp_df = pd.eval(f"{f}="+eq, target=decomp_df, engine='python')
            
        return decomp_df
        
    def _translate_string(self, s, d):
        return ''.join([d.get(x, x) for x in s])


#split apart for merging
df_1990 = df.loc[df["year_id"]==1990].reset_index()
df_2019 = df.loc[df["year_id"]==2019].reset_index()
df_2021 = df.loc[df["year_id"]==2021].reset_index()


#Component proportions
components = ["population", "sepsis_death_rate", 
              "amr_syns_prop", "bacterial_syns_prop", "res_bacterial_prop"]

das_gupta = das_gupta_decomp(N=5)

df_1990_2019 = das_gupta.decompose(base_df=df_1990, comp_df=df_2019, col_names = components)

component_effects = ["population_effect",
                    "sepsis_death_rate_effect",
                    "amr_syns_prop_effect",
                    "bacterial_syns_prop_effect",
                    "res_bacterial_prop_effect"]

df_1990_2019["calc_total_effect"] = df_1990_2019["population_effect"]+ df_1990_2019["sepsis_death_rate_effect"] + df_1990_2019["amr_syns_prop_effect"] + df_1990_2019["bacterial_syns_prop_effect"] + df_1990_2019["res_bacterial_prop_effect"]

merged_base_df = df_1990.merge(df_1990_2019, left_index=True, right_index=True, how='left')

grouped_df = merged_base_df[["location_id"] + components + component_effects].groupby(by=["location_id"]).agg({"population_effect":"sum",
                                                                         "sepsis_death_rate_effect": "sum",
                                                                         "amr_syns_prop_effect":"sum",
                                                                         "bacterial_syns_prop_effect":"sum",
                                                                         "res_bacterial_prop_effect":"sum"}).reset_index()

grouped_df['calc_tot_effect'] = grouped_df["population_effect"] + grouped_df["sepsis_death_rate_effect"] + grouped_df["amr_syns_prop_effect"] + grouped_df["bacterial_syns_prop_effect"] + grouped_df["res_bacterial_prop_effect"]

ab_df_1990 = assoc_burden.loc[assoc_burden["year_id"]==1990]
ab_df_2019 = assoc_burden.loc[assoc_burden["year_id"]==2019]
ab_df_2021 = assoc_burden.loc[assoc_burden["year_id"]==2021]

ab_df_1990_2019 = ab_df_1990.merge(ab_df_2019[["location_id", "amr_mean"]], how = "outer", on = "location_id")
ab_df_1990_2019["ab_diff"] = ab_df_1990_2019["amr_mean_y"]-ab_df_1990_2019["amr_mean_x"]

grouped_df = grouped_df.merge(ab_df_1990_2019[["location_id", "ab_diff"]], how = "left", on = "location_id")
df = grouped_df.copy()

for eff in component_effects:
    df[f"{eff}_perc"] = df[eff]/df["calc_tot_effect"] * 100
    
df["calc_perc"] = df["population_effect_perc"] + df["sepsis_death_rate_effect_perc"] + df["amr_syns_prop_effect_perc"] + df["bacterial_syns_prop_effect_perc"] + df["res_bacterial_prop_effect_perc"]

# To decompose (A) population change into (a1) population growth and (a2) change in population age structure, 
#we run a subsetted decomposition where we estimate the effect of A but rather than use 2019 populations each age group. 
#We then use the 1990 population * the all-age population growth from 1990 to 2019. 
#The difference between A and a1 will be a2

pop_df = pop_df.loc[pop_df["age_group_id"]==22]
pop_1990 = pop_df.loc[pop_df["year_id"]==1990]
pop_2019 = pop_df.loc[pop_df["year_id"]==2019]
pop_2021 = pop_df.loc[pop_df["year_id"]==2021]

#Prepping pop_growth df which has 1990, 2019, 2021
pop_growth = pop_1990.merge(pop_2019[["location_id", "population"]], how = "outer", on = "location_id")
pop_growth.rename(columns = {"population_x": "pop_1990",
                            "population_y": "pop_2019"}, inplace=True)
pop_growth = pop_growth.merge(pop_2021[["location_id", "population"]], how = "outer", on = "location_id")
pop_growth.rename(columns = {"population": "pop_2021"},inplace=True)
pop_growth = pop_growth[["location_id", "pop_1990", "pop_2019", "pop_2021"]]

pop_growth["pop_growth_rate_1990_2019"]= ((pop_growth["pop_2019"] - pop_growth["pop_1990"])/ pop_growth["pop_1990"]) + 1 
pop_growth["pop_growth_rate_2019_2021"]= ((pop_growth["pop_2021"] - pop_growth["pop_2019"])/ pop_growth["pop_2019"]) + 1 

#Decompose the population effect into population growth and age structure
pop_decomp = merged_base_df[["age_group_id", "location_id", "year_id"] + components]

pop_decomp.rename(columns={"population":"population_1990"}, inplace=True)
pop_decomp = pop_decomp.merge(df_2019[["age_group_id", "location_id", "population"]], how = "outer", on = ["age_group_id", "location_id"])
pop_decomp.rename(columns = {"population":"population_2019"}, inplace=True)

pop_decomp = pop_decomp.merge(pop_growth[["location_id", "pop_growth_rate_1990_2019"]], how = "left", on = "location_id")

#Apply the all age growth rate to each age group to get an alternative population for 2019
pop_decomp["all_age_rate_population_2019"] = pop_decomp["population_1990"] * pop_decomp["pop_growth_rate_1990_2019"]

#Prepare the alternative das gupta with the alternative population
pop_df_1990=pop_decomp.copy()
pop_df_1990["population"] = pop_df_1990["population_1990"]

#Replace population here with the alternative population (assuming 2019 population from the all age growth rate)
pop_df_2019 = df_2019.merge(pop_decomp[["location_id", "age_group_id", "all_age_rate_population_2019"]],
                           how = "left", on = ["location_id", "age_group_id"])

pop_df_2019.drop(columns=["population"], inplace=True)
pop_df_2019.rename(columns={"all_age_rate_population_2019":"population"}, inplace=True)

#Create the decomposition with the alternative 2019 population
df_1990_2019_pop_decomp = das_gupta.decompose(base_df=pop_df_1990, comp_df=pop_df_2019, col_names = components)

pop_merged_base_df = pop_df_1990.merge(df_1990_2019_pop_decomp, left_index=True, right_index=True, how='left')
pop_grouped_df = pop_merged_base_df[["location_id"] + components + component_effects].groupby(by=["location_id"]).agg({"population_effect":"sum",
                                                                         "sepsis_death_rate_effect": "sum",
                                                                         "amr_syns_prop_effect":"sum",
                                                                         "bacterial_syns_prop_effect":"sum",
                                                                         "res_bacterial_prop_effect":"sum"}).reset_index()
pop_grouped_df['calc_tot_effect'] = pop_grouped_df["population_effect"] + pop_grouped_df["sepsis_death_rate_effect"] + pop_grouped_df["amr_syns_prop_effect"] + pop_grouped_df["bacterial_syns_prop_effect"] + pop_grouped_df["res_bacterial_prop_effect"]

#Merge the normal full population effect with the counterfactual where the 2019 pop was based on the growth rate
comp_a_df = merged_base_df[["location_id", "age_group_id", "population_effect"]]
comp_a_df= comp_a_df.merge(pop_merged_base_df[["location_id", "age_group_id", "population_effect"]], how = "left", on = ["location_id", "age_group_id"])
comp_a_df.rename(columns={"population_effect_x":"pop_effect",
                         "population_effect_y": "population_growth_effect"}, inplace=True)

#Get the age structure change by taking the residual between the actual population effect and the population growth effect
comp_a_df["age_structure_change_effect"] = comp_a_df["pop_effect"] - comp_a_df["population_growth_effect"]
comp_a_df_grouped = comp_a_df[["location_id", "pop_effect", "population_growth_effect", "age_structure_change_effect"]].groupby(by=["location_id"]).agg({"pop_effect":"sum",
                                                                                                                                        "population_growth_effect":"sum",
                                                                                                                          "age_structure_change_effect":"sum"}).reset_index()

final_df = df.merge(comp_a_df_grouped, how ="left", on = "location_id")

component_effects_pop_decomp = [x for x in component_effects if x != "population_effect"]
component_effects_pop_decomp +=["population_growth_effect", "age_structure_change_effect"] 
final_df["test_total_effect"] = final_df[component_effects_pop_decomp].sum(axis=1)

final_df["population_growth_effect_perc"] = (final_df["population_growth_effect"]/final_df["calc_tot_effect"]) *100
final_df["age_structure_change_perc"] = (final_df["age_structure_change_effect"]/final_df["calc_tot_effect"])*100

final_df = add_location_metadata(final_df, ['location_name'])

full_decomp_1990_2019 = final_df[["location_id", "location_name", "population_effect", "population_growth_effect", 
                                  "age_structure_change_effect", "sepsis_death_rate_effect", 
                                  "amr_syns_prop_effect", "bacterial_syns_prop_effect", 
                                  "res_bacterial_prop_effect", "calc_tot_effect", "population_effect_perc", 
                                  "ab_diff", "population_growth_effect_perc", "age_structure_change_perc",
                                  "sepsis_death_rate_effect_perc", "amr_syns_prop_effect_perc", 
                                  "bacterial_syns_prop_effect_perc", "res_bacterial_prop_effect_perc", "calc_perc"]]

#full_decomp_1990_2019.to_parquet({FILEPATH})