% RUN_END2END_DEMO End-to-end demo for LCFV pipeline: SIFT -> GMM -> FV -> LCFV
%
%   Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
%   Signal Analysis and Machine Perception Laboratory,
%   Department of Electrical, Computer, and Systems Engineering,
%   Rensselaer Polytechnic Institute, Troy, NY 12180, USA
%
%   You are free to use this software for academic purposes if you cite our paper:
%   Quan Wang, Xin Shen, Meng Wang, Kim L. Boyer,
%   Label Consistent Fisher Vectors for Supervised Feature Aggregation,
%   22nd International Conference on Pattern Recognition (ICPR), 2014.
%
%   For commercial use, please contact the authors.

% Add paths
addpath(genpath(pwd));

clear; clc; close all;
fprintf('Running End-to-End LCFV Pipeline Demo\n');

%% 1. Generate Synthetic Data (Simulating SIFT descriptors)
% Assume 2 classes of images, 10 images each
% Each image has 100 descriptors of dim 64
n_classes = 2;
n_images_per_class = 10;
n_descs_per_image = 100;
D = 64;

fprintf('1. Generating synthetic SIFT-like descriptors...\n');
% Class 1: Centered at 0
% Class 2: Centered at 2
descs_pool = [];
labels = [];

% Train data for GMM (usually a separte large dataset, but reusing here for demo)
all_descs = [];
image_descs = {}; % Cell array to store descriptors for each image

img_idx = 1;
for c = 1:n_classes
    for i = 1:n_images_per_class
        % Generate descriptors
        base_mean = (c-1) * 2;
        descs = randn(D, n_descs_per_image) + base_mean;
        
        all_descs = [all_descs, descs];
        image_descs{img_idx} = descs;
        labels(img_idx) = c;
        img_idx = img_idx + 1;
    end
end

%% 2. Train GMM and PCA
K = 4; % Number of Gaussians (Small for demo speed)
pca_dim = 32; % Reduce dimension

fprintf('2. Training GMM (K=%d) and PCA (%d -> %d)...\n', K, D, pca_dim);
tic;
[w, mu, sigma, pca_transform, pca_mean] = fv_train(all_descs, K, pca_dim);
toc;

%% 3. Encode Fisher Vectors
fprintf('3. Encoding Fisher Vectors for %d images...\n', n_classes * n_images_per_class);
n_total_images = length(image_descs);
fv_dim = 2 * pca_dim * K;
FV = zeros(fv_dim, n_total_images);

tic;
for i = 1:n_total_images
    % Get raw descriptors
    X = image_descs{i};
    
    % Apply PCA
    X_centered = bsxfun(@minus, X, pca_mean);
    X_pca = pca_transform * X_centered;
    
    % Encode
    FV(:, i) = fv_encode(X_pca, w, mu, sigma);
end
toc;

fprintf('Fisher Vectors computed. Dimension: %d x %d\n', size(FV, 1), size(FV, 2));

%% 4. Apply LCFV
fprintf('4. Applying LCFV...\n');
% Ensure labels is a column vector
labels = labels(:);

% Create label comparison matrix
C1 = repmat(labels, 1, length(labels));
C2 = repmat(labels', length(labels), 1);
C = double(C1 == C2);

alpha = 10;
G = FV; % LCFV function expects D x N (same as our FV)

% LCFV1
tic;
[M1, W1] = solve_LCFV1(G, C, alpha);
LCFV1_Feats = M1 * G;
t1 = toc;
fprintf('LCFV1 computed in %.4f s\n', t1);

% LCFV2
tic;
M2 = solve_LCFV2(G, C, alpha);
LCFV2_Feats = M2 * G;
t2 = toc;
fprintf('LCFV2 computed in %.4f s\n', t2);

%% 5. Visualization (First 2 dimensions)
fprintf('5. Visualizing...\n');
f = figure('Visible', 'off'); % Invisible figure
% Helper to generate legend
legend_entries = arrayfun(@(c) sprintf('Class %d', c), 1:n_classes, 'UniformOutput', false);
colors = {'r', 'b', 'g', 'k', 'm', 'c'};
markers = {'o', 'x', '+', '*', 's', 'd'};

subplot(1, 3, 1);
hold on;
for c = 1:n_classes
    idx = (labels == c);
    plot(FV(1, idx), FV(2, idx), [colors{mod(c-1,6)+1} markers{mod(c-1,6)+1}]);
end
title('Original Fisher Vectors (First 2 Dim)');
xlabel('Dimension 1');
ylabel('Dimension 2');
legend(legend_entries, 'Location', 'best');
grid on;

subplot(1, 3, 2);
hold on;
for c = 1:n_classes
    idx = (labels == c);
    plot(LCFV1_Feats(1, idx), LCFV1_Feats(2, idx), [colors{mod(c-1,6)+1} markers{mod(c-1,6)+1}]);
end
title('LCFV1 (First 2 Dim)');
xlabel('Dimension 1');
ylabel('Dimension 2');
legend(legend_entries, 'Location', 'best');
grid on;

subplot(1, 3, 3);
hold on;
for c = 1:n_classes
    idx = (labels == c);
    plot(LCFV2_Feats(1, idx), LCFV2_Feats(2, idx), [colors{mod(c-1,6)+1} markers{mod(c-1,6)+1}]);
end
title('LCFV2 (First 2 Dim)');
xlabel('Dimension 1');
ylabel('Dimension 2');
legend(legend_entries, 'Location', 'best');
grid on;

fprintf('Saving visualization to synthetic_demo_results.png...\n');
print(f, 'synthetic_demo_results.png', '-dpng');
close(f);

fprintf('Demo completed successfully.\n');
