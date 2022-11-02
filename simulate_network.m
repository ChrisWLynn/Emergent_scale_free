% Script to simulate model of emergent scale-free networks, allowing self-
% loops and multi-edges. For efficiency, we only keep track of node
% degrees, not the entire adjacency matrix. We also compute the master
% equation at each step for comparison.

% Number of nodes:
N = 2*10^4;

% Number of edges:
E = 10^3;

% Average degree:
k_avg = 2*E/N;

% Probability of preferential attachment:
p = .5;

% Number of samples:
num_samples = 10^2 + 50;

% Number of steps per sample:
T_sample = N;

% Degrees to consider for master equation:
ks_master = (0:(2*E))';

%% Run calculation:

% Initialize with edges placed randomly:
G = generate_randNet(N, E);
ks = sum(G,2);

% Initialize distribution for master equation:
Pk_master = poisspdf(ks_master, k_avg);

% Quantities to record:
degrees = zeros(N, num_samples);
n = zeros(1, num_samples);
heterogeneity = zeros(1, num_samples);

% Master equation quantities to record:
deg_dist_master = zeros(sum(ks_master > 0), num_samples);
n_master = zeros(1, num_samples);
heterogeneity_master = zeros(1, num_samples);

% Record initial quantities:
degrees(:,1) = ks;
n(1) = sum(ks > 0);
heterogeneity(1) = 1/2*mean(abs(ks(ks > 0) - ks(ks > 0)'), [1 2])/mean(ks(ks > 0));

deg_dist_master(:,1) = Pk_master(ks_master > 0)/sum(Pk_master(ks_master > 0));
n_master(1) = N*(1 - Pk_master(1));
heterogeneity_master(1) = 1/2*sum((Pk_master(ks_master > 0)*Pk_master(ks_master > 0)')...
    .*abs(ks_master(ks_master > 0) - ks_master(ks_master > 0)'), [1 2])...
    /(sum(Pk_master(ks_master > 0))*sum(Pk_master(ks_master > 0).*ks_master(ks_master > 0)));

% Loop over samples:
for i = 2:num_samples
    
    tic
    
    % Loop over network updates per sample:
    counter = 1;
    
    while counter <= T_sample
    
        % Pick node to remove:
        I_remove = randi(N);
        while ks(I_remove) == 0
            I_remove = randi(N);
            counter = counter + 1;
            
            % Update master equation:
            Pk_master(2:end) = Pk_master(2:end) + 1/N*(-Pk_master(2:end)...
                + p*(-ks_master(2:end).*Pk_master(2:end) + ks_master(1:(end-1)).*Pk_master(1:(end-1)))...
                + (1-p)*k_avg*(-Pk_master(2:end) + Pk_master(1:(end-1))));
            Pk_master(1) = Pk_master(1) + 1/N*((-Pk_master(1) + 1) - (1-p)*k_avg*Pk_master(1));
            
        end
        
        % Remove edges:
        E_tot = ks(I_remove);
        ks(I_remove) = 0;
        
        % Number of preferential attachment updates:
        E_PA = binornd(E_tot, p);
        
        % Select nodes to connect via preferential attchment:
        Is_add_PA = randsample(N, E_PA, true, ks + realmin);
        
        % Select nodes to connect randomly:
        Is_add_rand = randsample(N, E_tot - E_PA, true);
        
        % Add new connections:
        Is_add = [Is_add_PA; Is_add_rand];
        for j = 1:E_tot
            ks(Is_add(j)) = ks(Is_add(j)) + 1;
        end
        
        % Increment counter:
        counter = counter + 1;
        
        % Update master equation:
        Pk_master(2:end) = Pk_master(2:end) + 1/N*(-Pk_master(2:end)...
            + p*(-ks_master(2:end).*Pk_master(2:end) + ks_master(1:(end-1)).*Pk_master(1:(end-1)))...
            + (1-p)*k_avg*(-Pk_master(2:end) + Pk_master(1:(end-1))));
        Pk_master(1) = Pk_master(1) + 1/N*((-Pk_master(1) + 1) - (1-p)*k_avg*Pk_master(1));
        
    end
    
    % Record quantities:
    degrees(:,i) = ks;
    n(i) = sum(ks > 0);
    heterogeneity(i) = 1/2*mean(abs(ks(ks > 0) - ks(ks > 0)'), [1 2])/mean(ks(ks > 0));
    
    % Record master equation quantities:
    deg_dist_master(:,i) = Pk_master(ks_master > 0)/sum(Pk_master(ks_master > 0));
    n_master(i) = N*(1 - Pk_master(1));
    heterogeneity_master(i) = 1/2*sum((Pk_master(ks_master > 0)*Pk_master(ks_master > 0)')...
        .*abs(ks_master(ks_master > 0) - ks_master(ks_master > 0)'), [1 2])...
        /(sum(Pk_master(ks_master > 0))*sum(Pk_master(ks_master > 0).*ks_master(ks_master > 0)));

end