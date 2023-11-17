%% Part 1
% Parameters
N = 1000;
p = 0.5;
q = 0.2;

% Generate SBM graph
dim = 2;
seed = 99;
rng(seed)
sig = randi([0, dim - 1], [N, 1]);
sig = 2 * rescale(sig) - 1;
P = (sig == sig') * p + (sig ~= sig') * q;
A = ones(N) - (rand(N) > P);
A = triu(A);
A = abs(A - A');

% Create graph and compute Fourier basis
G = gsp_graph(A);
G = gsp_create_laplacian(G, 'normalized');
G.type = 'normalized';
G = gsp_compute_fourier_basis(G);

% Choose subset of Fourier basis for visualization
U23 = G.U;
G.coords = U23(:, 2:max(3, dim));

% Plotting parameters
param = struct;
param.colorbar = 0;
G.plotting.edge_width = 0.1;
figure;
gsp_plot_signal(G, sig, param);
title('SBM Graph Clustering');
%% Part 2
% Parameters
N = 1000;
seed = 99;
dim = 2;

% Sets of alpha and beta
criteria = 0:0.05:1.5;
alphas = 4 * ones(1, numel(criteria));
betas = (sqrt(alphas) - sqrt(2*criteria)).^2;

% Initialize arrays to store results
acc = zeros(1, numel(criteria));

disp('Simulation start')

% Loop through sets
for i = 1:numel(alphas)
    alpha = alphas(i);
    beta = betas(i);
    
    % Generate SBM graph
    p = alpha * log(N) / N;
    q = beta * log(N) / N;

    % Run clustering and store accuracy
    acc(i) = run_clustering(N, p, q, dim, seed);
end

disp('Simulation finished')

% Find the point where accuracy becomes 1 after a certain criterion
criterion_threshold = 1;
idx_threshold = find(acc == 1, 1, 'first');

% Plot accuracy vs. criteria
figure;
plot(criteria, acc, '-o', 'LineWidth', 2);
hold on;
plot(criteria(idx_threshold), 1, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlabel('Criteria (\(\frac{(\sqrt{\alpha} - \sqrt{\beta})^2}{2}\))', 'Interpreter', 'latex');
ylabel('Accuracy');
title(['Accuracy vs. Criteria, N = ' num2str(N)]);
grid on;

% Display threshold information
disp(['Criterion threshold for N = ' num2str(N) ': ' num2str(criteria(idx_threshold))]);

function acc = run_clustering(N, p, q, dim, seed)
    % Generate SBM graph
    rng(seed)
    sig = randi([0, dim - 1], [N, 1]);
    sig = 2 * rescale(sig) - 1;
    P = (sig == sig') * p + (sig ~= sig') * q;
    A = ones(N) - (rand(N) > P);
    A = triu(A);
    A = abs(A - A');

    % Create graph and compute Fourier basis
    G = gsp_graph(A);
    G = gsp_create_laplacian(G, 'normalized');
    G.type = 'normalized';
    G = gsp_compute_fourier_basis(G);

    % Choose subset of Fourier basis for visualization
    U23 = G.U;
    dim = 2;
    G.coords = U23(:, 2:max(dim, 3));

    % Perform k-means clustering
    k = 2;
    [clusters, ~] = kmeans(G.coords, k);
    % clusters = G.coords(:,1) > 0;

    % True labels based on the sign of sig
    true_labels = (sig + 1) / 2 + 1;

    % Calculate accuracy
    accuracy = sum(clusters == true_labels) / numel(true_labels);
    accuracy = max(accuracy, 1 - accuracy);

    % Display results
    alpha = p * N / log(N);
    beta = q * N / log(N);
    c = (sqrt(alpha) - sqrt(beta)) ^ 2 / 2;
    disp(['alpha = ' num2str(alpha) ', beta = ' num2str(beta) ...
        ', criteria = ' num2str(c) ', Accuracy: ' num2str(accuracy)]);

    acc = accuracy;
    
end