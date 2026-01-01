% RUN_CIFAR10_END2END_DEMO Demo on a subset of CIFAR-10
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

addpath(genpath(pwd));
clear; clc; close all;

fprintf('Running CIFAR-10 Subset End-to-End Demo\n');
fprintf('This demo processes a subset of CIFAR-10 images using Dense SIFT + Fisher Vectors + LCFV.\n');

%% 1. Load Data
data_path = '../data/cifar-10-data_batch_1.mat';
if ~exist(data_path, 'file')
    error('CIFAR-10 data not found at %s', data_path);
end

fprintf('1. Loading Data from %s...\n', data_path);
load(data_path);
% 'data' is 10000x3072, 'labels' is 10000x1

% Settings
n_train = 50; % Small number for demo speed
n_test = 0;   % Just run on train for visualization demo
total_images = n_train;

% Select subset
raw_imgs = data(1:total_images, :);
labels_subset = labels(1:total_images);
% Fix labels to be 1-indexed for Matlab logic if needed, but CIFAR 0-9 is fine for logic
% Just ensuring labels is column
labels_subset = double(labels_subset(:)) + 1; % 1-10

%% 2. Extract Features (Dense SIFT)
fprintf('2. Extracting Dense SIFT Features from %d images...\n', total_images);
% Step size 4, Bin size 4 -> 16x16 patch
step_size = 8; % Larger step for speed in demo
bin_size = 4;

all_descs = [];
image_descs = {};

tic;
for i = 1:total_images
    % Reshape: 32x32x3
    % CIFAR format: row-major per channel? No, usually 1024 R, 1024 G, 1024 B
    img_flat = raw_imgs(i, :);
    R = reshape(img_flat(1:1024), 32, 32)';
    G = reshape(img_flat(1025:2048), 32, 32)';
    B = reshape(img_flat(2049:3072), 32, 32)';
    img = cat(3, R, G, B);
    img = uint8(img);
    
    descs = compute_dense_sift(img, step_size, bin_size);
    
    image_descs{i} = descs;
    all_descs = [all_descs, descs];
    
    if mod(i, 10) == 0
        fprintf('Processing image %d/%d...\n', i, total_images);
    end
end
t_extract = toc;
fprintf('Feature extraction took %.2f seconds\n', t_extract);
fprintf('Total descriptors: %d\n', size(all_descs, 2));

%% 3. Train GMM and PCA
K = 8; % Number of Gaussians
pca_dim = 64; 

fprintf('3. Training GMM (K=%d) and PCA (128 -> %d)...\n', K, pca_dim);
tic;
[w, mu, sigma, pca_transform, pca_mean] = fv_train(all_descs, K, pca_dim);
toc;

%% 4. Encode Fisher Vectors
fprintf('4. Encoding Fisher Vectors...\n');
fv_dim = 2 * pca_dim * K;
FV = zeros(fv_dim, total_images);

tic;
for i = 1:total_images
    X = image_descs{i};
    if isempty(X)
        % Handle case where no descriptors found (unlikely for dense)
        continue;
    end
    X_centered = bsxfun(@minus, X, pca_mean);
    X_pca = pca_transform * X_centered;
    FV(:, i) = fv_encode(X_pca, w, mu, sigma);
end
toc;

%% 5. Apply LCFV
fprintf('5. Applying LCFV...\n');
C1 = repmat(labels_subset, 1, length(labels_subset));
C2 = repmat(labels_subset', length(labels_subset), 1);
C = double(C1 == C2);

alpha = 10; % High alpha enforces stronger supervision

% LCFV1
[M1, W1] = solve_LCFV1(FV, C, alpha);
LCFV1_Feats = M1 * FV;

% LCFV2
M2 = solve_LCFV2(FV, C, alpha);
LCFV2_Feats = M2 * FV;

%% 6. Visualization
fprintf('6. Visualizing...\n');
f = figure('Name', 'CIFAR-10 Demo', 'Visible', 'off');

% Pick 5 classes to visualize to avoid clutter
classes_to_plot = unique(labels_subset);
if length(classes_to_plot) > 5
    classes_to_plot = classes_to_plot(1:5);
end

colors = {'r', 'b', 'g', 'k', 'm', 'c', 'y'};
markers = {'o', 'x', '+', '*', 's', 'd', '^'};

% Generate legend entries
legend_entries = {};
for k = 1:length(classes_to_plot)
    legend_entries{end+1} = sprintf('Class %d', classes_to_plot(k));
end

subplot(1, 3, 1);
hold on;
for k = 1:length(classes_to_plot)
    c = classes_to_plot(k);
    idx = (labels_subset == c);
    plot(FV(1, idx), FV(2, idx), [colors{mod(k-1,7)+1} markers{mod(k-1,7)+1}]);
end
title('Fisher Vectors (Input)');
xlabel('Dimension 1');
ylabel('Dimension 2');
legend(legend_entries, 'Location', 'best');
grid on;

subplot(1, 3, 2);
hold on;
for k = 1:length(classes_to_plot)
    c = classes_to_plot(k);
    idx = (labels_subset == c);
    plot(LCFV1_Feats(1, idx), LCFV1_Feats(2, idx), [colors{mod(k-1,7)+1} markers{mod(k-1,7)+1}]);
end
title('LCFV1 (Output)');
xlabel('Dimension 1');
ylabel('Dimension 2');
legend(legend_entries, 'Location', 'best');
grid on;

subplot(1, 3, 3);
hold on;
for k = 1:length(classes_to_plot)
    c = classes_to_plot(k);
    idx = (labels_subset == c);
    plot(LCFV2_Feats(1, idx), LCFV2_Feats(2, idx), [colors{mod(k-1,7)+1} markers{mod(k-1,7)+1}]);
end
title('LCFV2 (Output)');
xlabel('Dimension 1');
ylabel('Dimension 2');
legend(legend_entries, 'Location', 'best');
grid on;

fprintf('Saving visualization to cifar10_demo_results.png...\n');
print(f, 'cifar10_demo_results.png', '-dpng');
close(f);

fprintf('CIFAR-10 Demo Completed.\n');
