% Script to analyze real networks after dividing them into snapshots with a
% given number of edges. We compute degree distributions, network sizes,
% heterogeneities, and fit model parameters

% Load network information:
net_info = readtable('network_info.csv');

% Index of network to analyze:
net_ind = 1;

% Number of edges in each network snapshot:
E = 10^3;

% Load network:
edges = readmatrix(['Networks/', net_info.name{net_ind}, '.txt']);
    
% Number of edges:
num_edges = size(edges,1);

% Number of network snapshots:
num_snaps = floor(num_edges/E);

% List of snapshot labels:
net_inds = ceil((1:num_edges)/E);

% Degrees to consider:
ks = 1:(2*E);

% Record degree counts:
k_counts = zeros(1, length(ks));

% Other things to compute:
n = zeros(1, num_snaps);
heterogeneity = zeros(1, num_snaps);

% Loop over network snapshots:
for i = 1:num_snaps

    % Nodes in this snapshot:
    nodes = [edges(net_inds == i, 1); edges(net_inds == i, 2)];

    % Unique nodes:
    [nodes_unique, ~, nodes_new] = unique(nodes);
    n(i) = length(nodes_unique);

    % Degrees for this snapshot:
    ks_temp = histcounts(nodes, [nodes_unique', inf]);

    % Update counts:
    k_counts = k_counts + histcounts(ks_temp, [ks, inf]);

    % Compute degree heterogeneity:
    heterogeneity(i) = 1/2*(sum(abs(ks_temp' - ks_temp), [1 2])/...
        (n(i)*(n(i) - 1)))/mean(ks_temp);

end

% Compute degree distribution:
Pk = k_counts/sum(k_counts);

% Perform one-parameter model fit for preferential attachment probability p
% (NOTE: May need to adjust step size in fit_p):
k_avg = 2*E/net_info.num_nodes(net_ind);
p = fit_p(ks(Pk > 0), Pk(Pk > 0), k_avg, 0.5, 10^5);
