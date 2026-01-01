function [w, mu, sigma] = gmm_em(X, K, max_iters)
%GMM_EM Train Gaussian Mixture Model using Expectation-Maximization
%
%   [w, mu, sigma] = GMM_EM(X, K, max_iters)
%
%   Inputs:
%       X         - D x N data matrix (D dimensions, N samples)
%       K         - Number of Gaussians
%       max_iters - Maximum number of iterations (default: 100)
%
%   Outputs:
%       w         - 1 x K priors (weights)
%       mu        - D x K means
%       sigma     - D x K diagonal covariances (variances)
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

    if nargin < 3
        max_iters = 100;
    end

    [D, N] = size(X);
    
    % Initialization (Random subset of data)
    rand_idx = randperm(N, K);
    mu = X(:, rand_idx);
    
    % Initialize variances to be global variance
    global_var = var(X, 0, 2);
    sigma = repmat(global_var, 1, K);
    
    w = ones(1, K) / K;
    
    min_sigma = 1e-6; % Regularization
    
    last_ll = -inf;
    
    for iter = 1:max_iters
        % E-step: Compute responsibilities
        % P(z=k|x) propto w_k * N(x|mu_k, sigma_k)
        
        % Compute log likelihoods to avoid underflow
        log_rho = zeros(K, N);
        for k = 1:K
            % Log Gaussian: -0.5*D*log(2pi) - 0.5*sum(log(sigma_k)) - 0.5*sum((x-mu_k).^2 ./ sigma_k)
            diff = bsxfun(@minus, X, mu(:, k));
            sq_diff_norm = sum((diff.^2) ./ sigma(:, k), 1);
            log_det = sum(log(sigma(:, k)));
            log_rho(k, :) = log(w(k)) - 0.5 * (D * log(2 * pi) + log_det + sq_diff_norm);
        end
        
        % Normalize responsibilities using log-sum-exp trick
        max_log_rho = max(log_rho, [], 1);
        rho = exp(bsxfun(@minus, log_rho, max_log_rho));
        sum_rho = sum(rho, 1);
        gamma = bsxfun(@rdivide, rho, sum_rho);
        
        current_ll = sum(log(sum_rho) + max_log_rho);
        
        if iter > 1 && abs(current_ll - last_ll) < 1e-4 * abs(last_ll)
            fprintf('GMM converged at iteration %d\n', iter);
            break;
        end
        last_ll = current_ll;
        
        % M-step: Update parameters
        Nk = sum(gamma, 2)'; % 1 x K
        
        w = Nk / N;
        
        for k = 1:K
            % mu_k
            mu(:, k) = (X * gamma(k, :)') / Nk(k);
            
            % sigma_k
            diff = bsxfun(@minus, X, mu(:, k));
            sigma(:, k) = ((diff.^2) * gamma(k, :)') / Nk(k);
            
            % Regularize sigma
            sigma(:, k) = max(sigma(:, k), min_sigma);
        end
    end
end
