function fv = fv_encode(X, w, mu, sigma)
%FV_ENCODE Compute Fisher Vector for a set of descriptors
%
%   fv = FV_ENCODE(X, w, mu, sigma)
%
%   Inputs:
%       X      - D x N data matrix (descriptors)
%       w      - 1 x K priors
%       mu     - D x K means
%       sigma  - D x K diagonal covariances
%
%   Outputs:
%       fv     - (2*D*K) x 1 Fisher Vector (signed square rooted and L2 normalized)
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
    K = length(w);
    
    % Compute posteriors gamma(k, i)
    log_rho = zeros(K, N);
    
    % Precompute log terms
    log_w = log(w);
    
    for k = 1:K
        diff = bsxfun(@minus, X, mu(:, k));
        sq_diff_norm = sum((diff.^2) ./ sigma(:, k), 1);
        log_det = sum(log(sigma(:, k)));
        log_rho(k, :) = log_w(k) - 0.5 * (D * log(2 * pi) + log_det + sq_diff_norm);
    end
    
    max_log_rho = max(log_rho, [], 1);
    rho = exp(bsxfun(@minus, log_rho, max_log_rho));
    sum_rho = sum(rho, 1);
    gamma = bsxfun(@rdivide, rho, sum_rho);
    
    % In case of numerical issues
    gamma(isnan(gamma)) = 1/K;
    
    % Accumulate statistics
    % fv_u: gradients w.r.t mean
    % fv_v: gradients w.r.t variance
    
    fv_u = zeros(D, K);
    fv_v = zeros(D, K);
    
    for k = 1:K
        % D x N
        diff = bsxfun(@minus, X, mu(:, k));
        
        % gamma(k, :) is 1 x N
        gamma_k = gamma(k, :);
        
        S0 = sum(gamma_k);
        S1 = diff * gamma_k'; % D x 1
        S2 = (diff.^2) * gamma_k'; % D x 1
        
        fv_u(:, k) = (S1) ./ (sqrt(w(k)) * sqrt(sigma(:, k)));
        fv_v(:, k) = (S2 - S0 * sigma(:, k)) ./ (sqrt(w(k)) * sqrt(2) * sigma(:, k));
    end
    
    % Concatenate
    fv = [fv_u(:); fv_v(:)];
    
    % Power normalization
    fv = sign(fv) .* sqrt(abs(fv));
    
    % L2 normalization
    fv_norm = norm(fv);
    if fv_norm > 0
        fv = fv / fv_norm;
    end
end
