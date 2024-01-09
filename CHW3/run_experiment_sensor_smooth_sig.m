function [SNR_vals, SNR_avg] = run_experiment_sensor_smooth_sig(N, k)
    G = cell(k, 1);
    X = zeros(N, k);
    for i = 1:k
        G{i} = gsp_sensor(N);
        G{i} = gsp_compute_fourier_basis(G{i});
        alpha = rand(2, 1);
        X(:, i) = G{i}.U(:, 2:3) * alpha;
        scale = max(abs(max(X(:, i))), abs(min(X(:, i))));
        X(:, i) = X(:, i) / scale;
    end
    x = sum(X, 2);
    X_hat = optim_bss(G, x);
    SNR_vals = zeros(k, 1); 
    for i = 1:k
        SNR_vals(i) = snr(X(:, i), X_hat(:, i));
    end
    SNR_avg = mean(SNR_vals);
end