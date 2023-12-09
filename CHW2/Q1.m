%% Part 1: Create graph
seed = 140;
rng(seed)
% Define graph
W = [0 1 1 0 0 0 0 1;
      0 0 1 1 1 0 0 1;
      0 0 0 1 0 0 0 0;
      0 0 0 0 1 1 0 1;
      0 0 0 0 0 1 1 1;
      0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 1;
      0 0 0 0 0 0 0 0];
W = W + W';
param = struct();
param.distribute = 1;
G = gsp_sensor(8, param);
G.W = W;
G.plotting.vertex_size = 20;
% Plot graph
figure('Position', [100, 100, 800, 400]);
gsp_plot_graph(G);
title('Graph');
%% Part 2: Create original signal
G = gsp_create_laplacian(G);
G = gsp_compute_fourier_basis(G);
u1 = G.U(:, 1);
u2 = G.U(:, 2);
x = 2 * u1 + u2;
%% Part 3: Add Noise and Plot Signals
x_n = awgn(x, 10);
snr_orig = db(snr(x_n, x));
% Plotting
figure('Position', [100, 100, 1000, 600]);

subplot(2, 2, 1);
gsp_plot_signal(G, x);
title('Original Signal on Graph')

subplot(2, 2, 2);
gsp_plot_signal(G, x_n);
title('Noisy Signal on Graph', ['SNR: ', num2str(snr_orig, 4)])

% Plot the signals as time series
y_max = max([x; x_n]) + 0.1;
y_min = min(0, min([x; x_n]) - 0.1);

stem_subplot(x, 2, 2, 'Original Signal', 9, y_min, y_max, 3);
stem_subplot(x_n, 2, 2, 'Noisy Signal', 9, y_min, y_max, 4);

% Add an overall title for the subplots
sgtitle('Signal and Noisy Signal', 'FontSize', 16);

%% Part 4: Compute GFT, and Plot Spectra
Gn = G;
Gn = gsp_create_laplacian(Gn, 'normalized');
Gn.L = eye(Gn.N) - Gn.L;
Gn = gsp_compute_fourier_basis(Gn);

% Compute the GFT for all signals for both Laplacian and the new shift operator
x_hat_L = gsp_gft(G, x);
x_n_hat_L = gsp_gft(G, x_n);
x_hat_Wn = gsp_gft(Gn, x);
x_n_hat_Wn = gsp_gft(Gn, x_n);

% Plotting
figure('Position', [100, 100, 1000, 600]);

subplot(2, 2, 1);
gsp_plot_signal_spectral(G, x_hat_L);
title('Signal (Laplacian)')

subplot(2, 2, 2);
gsp_plot_signal_spectral(G, x_n_hat_L);
title('Noisy Signal (Laplacian)')

x_min = min(Gn.e);
x_max = max(Gn.e);

subplot(2, 2, 3);
gsp_plot_signal_spectral(Gn, x_hat_Wn);
xlim([x_min, x_max])
title('Signal (Norm-Weight)')

subplot(2, 2, 4);
gsp_plot_signal_spectral(Gn, x_n_hat_Wn);
xlim([x_min, x_max])
title('Noisy Signal (Norm-Weight)')

% Add an overall title for the subplots
sgtitle('Signal and Noisy Signal Spectrum', 'FontSize', 16);

%% Part 5: Creating ideal filters
filt_L =  [1 1 0 0 0 0 0 0]';
filt_Wn = [0 0 0 0 0 0 1 1]';

% Plotting
figure('Position', [100, 100, 1000, 400]);

subplot(1, 2, 1);
stem_filter(G.e, filt_L, 'Laplacian Filter');
xlabel('Frequency');
ylabel('Amplitude');

subplot(1, 2, 2);
stem_filter(Gn.e, filt_Wn, 'Norm-Weight Filter');
xlabel('Frequency');
ylabel('Amplitude');

%% Part 6: Filtering
x_filt_L = G.U * (filt_L .* x_n_hat_L);
x_filt_Wn = Gn.U * (filt_Wn .* x_n_hat_Wn);

% Calculate SNR
snr_filt_L = db(snr(x_filt_L, x));
snr_filt_Wn = db(snr(x_filt_Wn, x));

% Plotting
figure('Position', [100, 100, 1400, 600]);

subplot(2, 4, 1);
gsp_plot_signal(G, x);
title('Original Signal on Graph')

subplot(2, 4, 2);
gsp_plot_signal(G, x_n);
title('Noisy Signal on Graph', ['SNR: ', num2str(snr_orig, 4)])

subplot(2, 4, 3);
gsp_plot_signal(G, x_filt_L);
title('Filtered Signal(Laplacian) on Graph', ['SNR: ' num2str(snr_filt_L, 4)])

subplot(2, 4, 4);
gsp_plot_signal(G, x_filt_Wn);
title('Filtered Signal(Norm-Weight) on Graph', ['SNR: ' num2str(snr_filt_Wn, 4)])

% Plot the signals as time series
y_max = max([x; x_n; x_filt_Wn; x_filt_L]) + 0.1;
y_min = min(0, min([x; x_n; x_filt_Wn; x_filt_L]) - 0.1);

stem_subplot(x, 2, 4, 'Original Signal', 9, y_min, y_max, 5);
stem_subplot(x_n, 2, 4, 'Noisy Signal', 9, y_min, y_max, 6);
stem_subplot(x_filt_L, 2, 4, 'Filtered Signal(Laplacian)', 9, y_min, y_max, 7);
stem_subplot(x_filt_Wn, 2, 4, 'Filtered Signal(Norm-Weight)', 9, y_min, y_max, 8);

% Add an overall title for the subplots
sgtitle('Signal, Noisy Signal and Filtered Signals', 'FontSize', 16);
%% Part 7: Calculate SNR
% snr calculated in previous part. shown on graphs
%% Part 8: Best FIR fit
n = 3;
S123_L = G.e .^ (0:n-1);
S123_Wn =  Gn.e .^ (0:n-1);
h_L = pinv(S123_L) * filt_L;
h_L
FIR_L = S123_L * h_L;
h_Wn = pinv(S123_Wn) * filt_Wn;
h_Wn
FIR_Wn = S123_Wn * h_Wn;

% Plotting
figure('Position', [100, 100, 1000, 400]);

subplot(1, 2, 1);
stem(G.e, FIR_L, 'LineWidth', 1.5, 'DisplayName', 'FIR_L');
hold on;
stem(G.e, filt_L, 'LineWidth', 1.5, 'DisplayName', 'Ideal Filter');
hold off;
xlabel('Eigenvalue');
ylabel('Amplitude');
title('FIR Filter vs Ideal Filter (Laplacian)');
legend('Location', 'northeast');
grid on;

subplot(1, 2, 2);
stem(Gn.e, FIR_Wn, 'LineWidth', 1.5, 'DisplayName', 'FIR_{Wn}');
hold on;
stem(Gn.e, filt_Wn, 'LineWidth', 1.5, 'DisplayName', 'Ideal Filter');
hold off;
xlabel('Eigenvalue');
ylabel('Amplitude');
title('FIR Filter vs Ideal Filter (Norm-Weight)');
legend('Location', 'northwest');
grid on;

%% Part 9: Testing FIR filters
x_FIR_L = G.U * (FIR_L .* x_n_hat_L);
x_FIR_Wn = Gn.U * (FIR_Wn .* x_n_hat_Wn);

% Calculate SNR
snr_FIR_L = db(snr(x_FIR_L, x));
snr_FIR_Wn = db(snr(x_FIR_Wn, x));

% Plotting
figure('Position', [100, 100, 1400, 600]);

subplot(2, 4, 1);
gsp_plot_signal(G, x);
title('Original Signal on Graph')

subplot(2, 4, 2);
gsp_plot_signal(G, x_n);
title('Noisy Signal on Graph', ['SNR: ', num2str(snr_orig, 4)])

subplot(2, 4, 3);
gsp_plot_signal(G, x_FIR_L);
title('FIR-Filtered Signal(Laplacian) on Graph', ['SNR: ' num2str(snr_FIR_L, 4)])

subplot(2, 4, 4);
gsp_plot_signal(G, x_FIR_Wn);
title('FIR-Filtered Signal(Norm-Weight) on Graph', ['SNR: ' num2str(snr_FIR_Wn, 4)])

% Plot the signals as time series
y_max = max([x; x_n; x_FIR_Wn; x_FIR_L]) + 0.1;
y_min = min(0, min([x; x_n; x_FIR_Wn; x_FIR_L]) - 0.1);

stem_subplot(x, 2, 4, 'Original Signal', 9, y_min, y_max, 5);
stem_subplot(x_n, 2, 4, 'Noisy Signal', 9, y_min, y_max, 6);
stem_subplot(x_FIR_L, 2, 4, 'FIR-Filtered Signal(Laplacian)', 9, y_min, y_max, 7);
stem_subplot(x_FIR_Wn, 2, 4, 'FIR-Filtered Signal(Norm-Weight)', 9, y_min, y_max, 8);

% Add an overall title for the subplots
sgtitle('Signal, Noisy Signal and FIR-Filtered Signals', 'FontSize', 16);
%% Helper functions
% Function to create a stem plot in a subplot
function stem_subplot(x, n, m, title_str, xlim_val, y_min, y_max, subplot_pos)
    subplot(n, m, subplot_pos);
    stem(x);
    title(title_str);
    xlim([0, xlim_val]);
    ylim([y_min, y_max]);
end

function stem_filter(eigenvalues, filter, title_str)
    stem(eigenvalues, filter);
    hold on
%     plot(eigenvalues, filter);
    title(title_str);
end