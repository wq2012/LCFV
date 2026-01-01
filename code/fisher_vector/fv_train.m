function [w, mu, sigma, pca_transform, pca_mean] = fv_train(X, K, pca_dim)
%FV_TRAIN Train GMM and PCA for Fisher Vector encoding
%
%   [w, mu, sigma, pca_transform, pca_mean] = FV_TRAIN(X, K, pca_dim)
%
%   Inputs:
%       X             - D x N data matrix (descriptors)
%       K             - Number of Gaussians
%       pca_dim       - Target dimension for PCA (optional, if empty no PCA)
%
%   Outputs:
%       w, mu, sigma  - GMM parameters
%       pca_transform - PCA projection matrix (pca_dim x D)
%       pca_mean      - Mean descriptor (D x 1)
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

    [D, N] = size(X);
    
    % PCA implementation
    if nargin > 2 && ~isempty(pca_dim) && pca_dim < D
        fprintf('Training PCA (%d -> %d)...\n', D, pca_dim);
        pca_mean = mean(X, 2);
        X_centered = bsxfun(@minus, X, pca_mean);
        
        % SVD on covariance is better for memory than SVD on X if N >> D
        C = (X_centered * X_centered') / (N - 1);
        [U, S, ~] = svd(C);
        
        pca_transform = U(:, 1:pca_dim)';
        
        % Project data for GMM training
        X = pca_transform * X_centered;
    else
        pca_mean = zeros(D, 1);
        pca_transform = eye(D);
    end
    
    fprintf('Training GMM (K=%d)...\n', K);
    [w, mu, sigma] = gmm_em(X, K);
end
