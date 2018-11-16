%%  A modified version of MATLAB cholcov() function
%   Cholesky-like covariance decomposition

function [T,p] = cholcov2(Sigma)

[n,m] = size(Sigma);
if n~=m
    error('Sigma is not square');
end

Sigma=(Sigma+Sigma')/2;

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
        D(D<0)=0;
        T = U(:,t) * diag(sqrt(D)) * U(:,t)';
    end
end

