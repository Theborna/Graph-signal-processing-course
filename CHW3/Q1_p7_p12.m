%% Part 7: Loading images and formatting
face1 = imread('face1.jpg');
face2 = imread('face2.jpg');

% Transform images to double grayscale
face1 = double(rgb2gray(face1));
face2 = double(rgb2gray(face2));

% Vectorize images
x1 = reshape(face1, [], 1);
x2 = reshape(face2, [], 1);
X = [x1, x2];

%% Part 8: Combining signals
% Combine images
x = x1 + x2;

% Reshape to image dimensions
combined_image = reshape(x / 2, size(face1));

%% Visualization
figure('Position', [100, 100, 800, 300]);

subplot(1, 3, 1);
imshow(face1, []);
title('Face 1');

subplot(1, 3, 2);
imshow(combined_image, []);
title('Combined Image (x1 + x2)');

subplot(1, 3, 3);
imshow(face2, []);
title('Face 2');

% Visual divider
disp(repmat('=', 1, 50));

%% Part 9: Creating appropriate graphs
G = cell(2, 1);
[N, M] = size(face1);

for i = 1:2
    % Create Kings graph
    param.rule = 'strong';
    kings = gsp_graph_product(gsp_path(N), gsp_path(M), param);
    G{i} = gsp_2dgrid(N, M);
    G{i}.W = kings.W;

    % Modify graph weights based on image data
    dij = abs(X(:, i)' - X(:, i));
    eta = 0.001;
    wij = 1 ./ (dij + eta);
    G{i}.W = G{i}.W .* wij;
    
    G{i} = gsp_create_laplacian(G{i}, 'combinatorial');
    % Compute Fourier basis
    G{i} = gsp_compute_fourier_basis(G{i});
end
%% Part 10: Visualizing smoothness
figure('Position', [100, 100, 800, 400]);

t = 500;
for i = 1:2
    subplot(1, 2, i);
    
    % Perform GFT on the signal
    x_hat = gsp_gft(G{i}, X(:, i));
    
    % Plot GFT coefficients
    stem(G{i}.e(1:t), x_hat(1:t), 'Marker', 'none', 'LineWidth', 1.5);
    
    title(['Graph ' num2str(i)]);
    xlabel('Frequency');
    ylabel('Amplitude');
    xlim([G{i}.e(1), G{i}.e(t)]);
    grid on;
end

%% Part 11: BSS
x = x - mean(x);
X_hat = optim_bss(G, x);
x1_recon = X_hat(:, 1);
x2_recon = X_hat(:, 2);

%% Calculate SNR
snr_1 = snr(x1_recon, x1 - mean(x1));
snr_2 = snr(x2_recon, x2 - mean(x2));

%% Visualization
face1_recon = reshape(x1_recon, size(face1));
face2_recon = reshape(x2_recon, size(face2));
face1_recon = rescale(face1_recon, 0, 255);
face2_recon = rescale(face2_recon, 0, 255);
figure('Position', [100, 100, 1200, 400]);

% Plotting reconstructed faces
subplot(2, 3, 1);
imshow(face1_recon, []);
title(['Reconstructed Face 1 (SNR: ' num2str(snr_1, '%.2f') ' dB)']);

subplot(2, 3, [2, 5]);
imshow(combined_image, []);
title('Combined Image (x1 + x2)');

subplot(2, 3, 3);
imshow(face2_recon, []);
title(['Reconstructed Face 2 (SNR: ' num2str(snr_2, '%.2f') ' dB)']);

subplot(2, 3, 4);
imshow(face1, []);
title('Face 1');

subplot(2, 3, 6);
imshow(face2, []);
title('Face 2');

% Visual divider
disp(repmat('=', 1, 50));