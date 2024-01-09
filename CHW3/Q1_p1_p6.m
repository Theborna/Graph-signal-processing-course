%% Part 1 - Part 4
k = 4;
[SNR_i, SNR_avg] = run_experiment_sensor_smooth_sig(100, k);

%% Part 5
k = 4;
SNR_avg_i = zeros(1, 20);
for i = 1:20
    [~, SNR_avg_i(i)] = run_experiment_sensor_smooth_sig(100, k);
end
true_SNR_avg = mean(SNR_avg_i);

%% Part 6
k = 4;
SNR_avg_i_200 = zeros(1, 20);
SNR_avg_i_300 = zeros(1, 20);
for i = 1:20
    [~, SNR_avg_i_200(i)] = run_experiment_sensor_smooth_sig(200, k);
    [~, SNR_avg_i_300(i)] = run_experiment_sensor_smooth_sig(300, k);
end
true_SNR_avg_200 = mean(SNR_avg_i_200);
true_SNR_avg_300 = mean(SNR_avg_i_300);

%% Unimportant (Boxplot)
data = [SNR_avg_i', SNR_avg_i_200', SNR_avg_i_300'];
figure('Position',[100, 100, 800, 600]);

boxplot(data, 'Labels', {'N=100', 'N=200', 'N=300'}, 'Notch','on');
title('Comparison of Signal-to-Noise Ratios (SNR) at Different dimensionality');
xlabel('Smooth Signal dimention');
ylabel('Signal-to-Noise Ratio (SNR)');
grid on;

disp(repmat('=', 1, 50));
