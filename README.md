# Emergent_scale_free
Code used to perform the analyses in the paper "Emergent scale-free networks" by Christopher W. Lynn, Caroline M. Holmes, and Stephanie E. Palmer.

The real networks are listed in "network_info.csv". The data for each network can be found on Dropbox:

For each network, each row corresponds to a distinct edge. The edges are listed in the temporal order in which they occurred.

Here, we include the following scripts and functions:

"simulate_network.m": Script to stimulate our model and compute the evolution of the master equation.

"analyze_real_networks.m": Script to analyze the real networks.

"deg_dist_analytic.m": Function to compute the analytic degree distribution in our model.

"fit_p.m": Function used to fit the preferential attachment probability p in our model to a real network (used in "analyze_real_networks.m").

"generate_randNet.m": Function to generate a random network (used in "simulate_network.m").
