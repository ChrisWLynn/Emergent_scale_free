function [G] = generate_randNet(N, E)
% Inputs: number of nodes N and number of edges E
%
% Ouput: NxN undirected adjacency matrix G with E edges randomly placed
% among the nodes. We allow self-loops and multi-edges.

% Two random lists of nodes:
Is = randi(N, E, 1);
Js = randi(N, E, 1);

% Create matrix:
G = full(sparse(Is, Js, ones(E,1), N, N));
G = G + G';
