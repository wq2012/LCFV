function [T,p] = cholcov2(Sigma)
%CHOLCOV2 A modified version of MATLAB cholcov() function
%   Cholesky-like covariance decomposition
%
%   [T,p] = CHOLCOV2(SIGMA) computes the Cholesky-like decomposition of the
%   covariance matrix SIGMA.
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

    [n,m] = size(Sigma);
    if n~=m
        error('Sigma is not square');
    end

    Sigma = (Sigma + Sigma') / 2;

    [T,p] = chol(Sigma);

    if p > 0
        % Test for positive definiteness

        % Can get factors of the form Sigma==T'*T using the eigenvalue
        % decomposition of a symmetric matrix, so long as the matrix
        % is positive semi-definite.
        [U,D] = eig(full(Sigma));

        % Pick eigenvector direction so max abs coordinate is positive
        [~,maxind] = max(abs(U),[],1);
        negloc = (U(maxind + (0:n:(m-1)*n)) < 0);
        U(:,negloc) = -U(:,negloc);

        D = diag(D);
        tol = eps(max(D)) * length(D);
        t = (abs(D) > tol);
        D = D(t);
        p = sum(D<0); % number of negative eigenvalues

        if (p==0)
            T = U(:,t) * diag(sqrt(D)) * U(:,t)';
        else
            D(D<0) = 0;
            T = U(:,t) * diag(sqrt(D)) * U(:,t)';
        end
    end
end

