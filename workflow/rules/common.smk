import pandas as pd

# this file has one line per base simulation run
# up through the recap and mutate step
simulations = pd.read_csv(config["simulations"], sep="\t", comment='#')
assert len(simulations) == len(simulations['sim_name'].unique()), "simulation names (sim_name) must be unique."


# ascertainment
ascertainments = pd.read_csv(config["ascertainments"], sep="\t", comment='#')
assert len(ascertainments) == len(ascertainments['asc_name'].unique()), "ascertainment names (asc_name) must be unique."


# this file has one line per analysis run
# specifies the sampling and filtering
analyses = pd.read_csv(config["analyses"], sep="\t", comment='#')
assert len(analyses) == len(analyses['anal_name'].unique()), "analysis names (anal_name) must be unique."

units = simulations.merge(ascertainments, on='sim_name')

units = units.merge(analyses, on='asc_name')

# used for getting the names of the MOSAIC output files
units['nsource'] = units['target_pop']

simulations = simulations.set_index("sim_name", drop=False)
ascertainments = ascertainments.set_index("asc_name", drop=False)
analyses = analyses.set_index("anal_name", drop=False)
# merge them all - this is the dataframe used in most cases

units = units.set_index(['sim_name', 'asc_name', 'anal_name'], drop=False)
