function [p, Pk_pred, err, i] = fit_p(ks, Pk, k_avg, p0, max_steps)
% Inputs: Observed degrees ks, probabilities of different degrees Pk, and
% average degree k_avg. We also take an initial guess at the preferential
% attachment probability p0. We note that k_avg is the average strength
% over all possible nodes (including k = 0), and therefore will likely not
% equal sum(ks.*Pk).
%
% Outputs: Fit preferential attachment probability p to observed degree
% distribution. Also return predicted degree distribution Pk_pred, fit
% error err, and number of gradient descent steps i.

% Step size (NOTE: may need to adjust this for different networks):
step_size = 10^(-3);

% Derivative resolution:
dp = 10^(-10);

% Stopping criteria:
dE_dp_min = 10^(-6);

% Make sure we only have unique degrees:
ks_unique = unique(ks);
Pk_unique = zeros(size(ks_unique));
num_s = length(ks_unique);

for i = 1:num_s
    Pk_unique(i) = sum(Pk(ks == ks_unique(i)));
end

% Make sure probabilities are normalized to one:
Pk_unique = Pk_unique/sum(Pk_unique);

% Initialize preferential attachment probability:
p = p0;

% Loop until max number of steps:
for i = 1:max_steps
    
    % Compute model probabilities:
    Pk_pred = deg_dist_analytic(ks_unique, p, k_avg);
    
    % Compute model probabilities lightly away from p (for derivative):
    p_der = p + dp;
    Pk_der = deg_dist_analytic(ks_unique, p_der, k_avg);
    
    % Compute error and derivative of error:
    err = sqrt(sum((log(Pk_unique) - log(Pk_pred)).^2));
    err_der = sqrt(sum((log(Pk_unique) - log(Pk_der)).^2));
    
    % Compute derivatives of error:
    dE_dp = (err_der - err)/dp;
    
    % Stop if derivative is small enough:
    if abs(dE_dp) < dE_dp_min
        return;
    end
    
    % Update PA probability:
    p = p - step_size*dE_dp/sqrt(i);
    
    % Stop if p becomes ill-defined:
    if p < 0
        
        p = 0;
        Pk_pred = deg_dist_analytic(ks_unique, p, k_avg);
        err = sqrt(sum((log(Pk_unique) - log(Pk_pred)).^2));
        return;

    elseif p > 1

        p = 1;
        Pk_pred = deg_dist_analytic(ks_unique, p, k_avg);
        err = sqrt(sum((log(Pk_unique) - log(Pk_pred)).^2));
        return;

    end

end



