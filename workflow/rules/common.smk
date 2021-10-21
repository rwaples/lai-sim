import pandas as pd

# this file has one line per base simulation run
# up through the recap and mutate step
simulations = pd.read_csv(config["simulations"], sep="\t")
assert len(simulations == len(simulations['sim_name'].unique())), "simulation names (sim_name) must be unique."

# this file has one line per analysis run
# each line should reference a simulation
# specifies the sampling and filtering
# may also specify a limited genomic span?? (probably not worth it)
# we will have to see how best to specify the LAI parameters
analyses = pd.read_csv(config["analyses"], sep="\t")
assert len(analyses == len(analyses['anal_name'].unique())), "analysis names (anal_name) must be unique."

units = analyses.merge(simulations, on=['sim_name'])

# used for getting the names of the MOSAIC output files
units['nsource'] = units['target_pop']
units['naming_mess'] = units['nind_admixed'].map(str) + "_22-22_" + (2*(units['nind_admixed']+units['nsource']*units['nind_ref'])).map(str) + "_60_0.99_100"

simulations = simulations.set_index("sim_name", drop=False)
analyses = analyses.set_index("anal_name", drop=False)
units = units.set_index(['sim_name', 'anal_name'], drop = False)
