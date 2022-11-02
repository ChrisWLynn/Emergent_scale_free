function Pk = deg_dist_analytic(ks, p, k_avg)
% Inputs: list of degrees ks, preferential attachment probability p, and
% average degree k_avg
%
% Outputs: Analytic degree distribution Pk for our emergent scale-free
% model, evaluated at the degrees ks and normalized to run over all
% positive degrees k > 0.
%
% NOTE: The average degree k_avg represents the average over all degrees
% (including k = 0), which is different from sum(Pk.*ks), the average over
% positive degrees.

% Make sure we only consider positive degrees:
if min(ks) < 1
    error('Only consider positive degrees!');
end

% Check if p = 0:
if p == 0

    Pk = (k_avg.^(ks - 1))./((1 + k_avg).^ks);
    return;

end

% Compute probabilities:
Pk = gamma(1 + 1/p)*gamma(1 + 1/p + k_avg*(1/p - 1))/(gamma(1/p)*gamma(1 + k_avg*(1/p - 1)))...
    *gamma(ks + k_avg*(1/p - 1))./gamma(ks + 1 + 1/p + k_avg*(1/p - 1));

% Deal with degrees for which distribution is poorly defined:
inds_bad = union(find(isnan(Pk)), find(Pk < realmin));
if ~isempty(inds_bad)
    c = Pk(min(inds_bad)-1)/ks(min(inds_bad)-1)^(-(1+1/p));
    Pk(inds_bad) = c*ks(inds_bad).^(-(1+1/p));
end