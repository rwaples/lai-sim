import pandas as pd

# this file has one line per base simulation run
# up through the recap and mutate step
simulations = (pd.read_csv(config["simulations"], sep="\t")
    .set_index("sim_name", drop=False)
    .sort_index()
)

# this file has one line per analysis run
# each line should reference a simulation
# specifies the sampling and filtering
# may also specify a limited genomic span?? (probably not worth it)
# we will have to see how best to specify the LAI parameters
analyses = (pd.read_csv(config["analyses"], sep="\t")
    .set_index("name", drop=False)
    .sort_index()
)
