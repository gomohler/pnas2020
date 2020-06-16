# pnas2020

Branching process:

1.  Run “create_data_new_format.R” to get covid-19 case data from JHU GitHub.

2.  Run “run_state_confirmed.m” to train branching process on confirmed cases, run “run_states.m” to train branching process model on death counts.

3.  csv files “mod_compare_state” and “mod_compare_state_confirmed” contain data corresponding to table 1 in the paper.

Comparison to SIR/SEIR

1.  Run model_compare.R to create sir/seir data for table 1 in the paper.

Country level dynamic R estimate

1. Code located in dynamic_R folder (run create_data_updated.R followed by run_all_updated.m)

