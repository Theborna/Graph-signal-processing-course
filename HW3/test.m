N = 10;
M = 5;
G1 = gsp_full_connected(N);
G2 = gsp_full_connected(M);

H = gsp_graph_product(G1, G2);

H = gsp_compute_fourier_basis(H);

t = 0.3;

H.coords = t * kron(G1.coords, ones(M, 1)) + (1-t) * kron(ones(N,1), G2.coords);

% U = H.U;

U = kron(dftmtx(M)/sqrt(M), dftmtx(N)/sqrt(N));

MN = M * N;

Q = zeros(MN);
T = zeros(MN);

for l = 1:MN
    for m = 1:MN
        d_m = zeros(MN, 1);
        d_l = zeros(MN, 1);
        d_m(m) = 1;
        d_l(l) = 1;
        
        d_hat_m = (U') * d_m;
        d_hat_l = U' * d_l;
        
        y_hat = d_hat_m .* d_hat_l;
        
        y = U * y_hat;
        
        Q(l,m) = find(y == max(y));
        m1 = floor((m-1)/N);
        m2 = mod((m-1), N);
        l1 = floor((l-1)/N);
        l2 = mod(l-1, N);

        t2 = mod(m2+l2, N);
        t1 = mod(m1+l1, M);
        
        t0 = N*t1 + t2;

        T(l,m) = mod(t0, MN) + 1; 
    end
end

subplot(1,2,1)
surf(Q)
subplot(1,2,2)
surf(T)
% gsp_plot_signal(H, y)