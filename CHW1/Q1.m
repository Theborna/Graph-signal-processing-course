%% part 1
% Define graph 1
w1 = [0, 0.7, 1.1, 2.3];
W1 = zeros(4);
W1(1, :) = w1;
W1 = W1 + W1';
G1 = gsp_comet(4, 3);
G1.W = W1;
G1.plotting.limits = [-0.5, 1, -1, 1];
G1.plotting.vertex_size = 20;

% Define graph 2
W2 = [0, 1.6, 2.4; 0, 0, 0.8; 0, 0, 0];
W2 = W2 + W2';
G2 = gsp_ring(3);
G2.W = W2;
G2.coords = G2.coords - mean(G2.coords);
G2.plotting.limits = [-0.5, 1, -1, 1];
G2.plotting.vertex_size = 20;

% Plot graphs
figure('Position', [100, 100, 800, 400]);

subplot(1, 2, 1);
gsp_plot_graph(G1);
title('Graph 1');

subplot(1, 2, 2);
gsp_plot_graph(G2);
title('Graph 2');
%% part 2
pt.verbose = 1;
pt.rule = 'kronecker';
Gs = gsp_graph_product(G1, G2, pt);

Gs.coords
pt.rule = 'cartesian';
Gt = gsp_graph_product(G1, G2, pt);

alpha = 0.25;
Gs.coords = alpha * kron(ones(size(G1.coords, 1), 1), G2.coords) + ...
    (1 - alpha) * kron(G1.coords, ones(size(G2.coords, 1), 1));

Gt.coords = Gs.coords;

% Plot graphs
figure('Position', [100, 100, 800, 400]);

subplot(1, 2, 1);
gsp_plot_graph(Gs);
title('Kronecker Product Rule');

subplot(1, 2, 2);
gsp_plot_graph(Gt);
title('Cartesian Product Rule');

% Set a centered title for the entire figure
sgtitle('Graph Products: Kronecker vs. Cartesian', 'FontSize', 16);

% Display and save adjacency matrices
display_and_save_matrices(Gs, 'Gs_matrices.txt');
display_and_save_matrices(Gt, 'Gt_matrices.txt');

%% Part 3
H = Gt;
x = rand(size(H.A, 1), 1) * 20 - 10;

% Plot signal
figure('Position', [100, 100, 600, 500]);
title('Random signal on graph')
gsp_plot_signal(H, x)

%% Part 4
H  = gsp_compute_fourier_basis(H);
spectrum = H.e;
x_hat = gsp_gft(H, x);
figure('Position', [100, 100, 800, 400]);
% gsp_plot_signal_spectral(H,x_hat);
plot(real(spectrum), imag(spectrum), '*', 'MarkerSize', 5)
%% Part 5
v1 = H.U(:, 1);
v2 = H.U(:, 2);
v_end1 = H.U(:, end);
v_end2 = H.U(:, end-1);
% Plot signal
figure('Position', [100, 100, 900, 600]);

subplot(2, 2, 1);
gsp_plot_signal(H, v1);
title(['First Eigenvector, ' '$\lambda$ = ' num2str(round(H.e(1), 2))],...
    'Interpreter','latex');

subplot(2, 2, 2);
gsp_plot_signal(H, v2);
title(['Second Eigenvector, ' '$\lambda$ = ' num2str(round(H.e(2), 2))],...
    'Interpreter','latex');

subplot(2, 2, 3);
gsp_plot_signal(H, v_end1);
title(['Last Eigenvector, ' '$\lambda$ = ' num2str(round(H.e(end), 2))],...
    'Interpreter','latex');

subplot(2, 2, 4);
gsp_plot_signal(H, v_end2);
title(['One to the last Eigenvector, ' '$\lambda$ = ' num2str(round(H.e(end-1), 2))],...
    'Interpreter','latex');

% Add an overall title for the subplots
sgtitle('Graph Signal Analysis', 'FontSize', 16);
