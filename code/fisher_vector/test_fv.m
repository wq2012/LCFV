function test_fv()
%TEST_FV Unit tests for Fisher Vector module

    addpath(genpath(pwd));
    
    test_gmm_convergence();
    test_fv_encoding_shape();
    test_pca();
    
    fprintf('All FV tests passed.\n');
end

function test_gmm_convergence()
    fprintf('Testing GMM convergence...\n');
    % Synthetic data: 3 clusters
    N = 300;
    D = 2;
    mu_true = [0, 5, 10; 0, 5, 0];
    X = [bsxfun(@plus, randn(D, N/3), mu_true(:,1)), ...
         bsxfun(@plus, randn(D, N/3), mu_true(:,2)), ...
         bsxfun(@plus, randn(D, N/3), mu_true(:,3))];
     
    K = 3;
    [w, mu, sigma] = gmm_em(X, K, 50);
    
    assert(length(w) == K, 'Wrong weight dimension');
    assert(all(size(mu) == [D, K]), 'Wrong mean dimension');
    assert(all(size(sigma) == [D, K]), 'Wrong sigma dimension');
    % Roughly check if means are close to [0,0], [5,5], [10,0] (order may vary)
    % Just checking if code runs without error and returns valid values here
    assert(all(~isnan(mu(:))), 'Means contain NaN');
end

function test_fv_encoding_shape()
    fprintf('Testing FV encoding shape...\n');
    D = 10;
    K = 4;
    N = 50;
    
    X = rand(D, N);
    
    % Fake GMM
    w = ones(1, K)/K;
    mu = rand(D, K);
    sigma = ones(D, K);
    
    fv = fv_encode(X, w, mu, sigma);
    
    assert(all(size(fv) == [2*D*K, 1]), 'Wrong FV dimension');
    assert(abs(norm(fv) - 1) < 1e-5, 'FV should be L2 normalized');
end

function test_pca()
    fprintf('Testing PCA...\n');
    D = 10;
    D_new = 5;
    N = 100;
    X = randn(D, N);
    
    [~, ~, ~, pca_transform, ~] = fv_train(X, 1, D_new);
    
    assert(all(size(pca_transform) == [D_new, D]), 'Wrong PCA matrix dimension');
end
