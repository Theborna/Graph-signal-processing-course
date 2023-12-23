N = 8;
G1 = gsp_full_connected(N);
G2 = gsp_full_connected(N);

H = gsp_graph_product(G1, G2);

H = gsp_compute_fourier_basis(H);

t = 0.3;

H.coords = gsp_2dgrid(N, N).coords;

x = zeros(N*N, 1);
x(1) = 64;

x_hat = H.U' * x;
x_filt_hat = zeros(N*N, 1);
x_filt_hat(1:15) = x_hat(1:15);
x_filt = H.U * x_filt_hat;

x_proc = x;
for i = 1:20
    x_proc = H.L * x_proc;
    x_proc = x_proc / sqrt(x_proc'*x_proc);
end

figure('Position', [100, 100, 1200, 300])
subplot(1, 3, 1)
gsp_plot_signal(H, x)
title('original')
subplot(1, 3, 2)
gsp_plot_signal(H, x_filt)
title('filtered')
subplot(1, 3, 3)
gsp_plot_signal(H, x_proc)
title('process')