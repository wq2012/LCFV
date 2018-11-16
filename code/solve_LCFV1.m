% Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA
% 
% You are free to use this software for academic purposes if you cite our paper: 
% Quan Wang, Xin Shen, Meng Wang, Kim L. Boyer, 
% Label Consistent Fisher Vectors for Supervised Feature Aggregation, 
% 22nd International Conference on Pattern Recognition (ICPR), 2014. 
% 
% For commercial use, please contact the authors. 


%%  This is the implementation of LCFV1
%   Method: 
%       G'*W*G=G'*G+C, and W=I+E, then solve E in G'*E*G=C, 
%       then W=M'*M
%   G: m*n Fisher vectors
%   C: n*n label comparison matrix
%   alpha: parameter
%   M: m*m resulting transformation matrix to be applied on Fisher vectors
%   W: m*m matrix, W=M'*M
%   M*G is called LCFV
%   m: dimension of features
%   n: number of instances

function [M,W]=solve_LCFV1(G,C,alpha)

[m,n]=size(G);

C=C*alpha;

sv_t=0.01; % threshold for small singular values

if n>=m % over-determined
    disp('LCFV1: over-determined problem');
    E = ( G' \ C ) / G;
else % under-determined
    disp('LCFV1: under-determined problem');
    [U,S,V] = svd(G); % X = U*S*V'
    nn=sum(diag(S)>sv_t); % remove small singular values
    S=S(:,1:nn); % remove small singular values
    V=V(:,1:nn); % remove small singular values
    C=V'*C*V;
    Sinv=diag(1./diag(S));
    Z1=Sinv*C*Sinv;
    Z=[Z1 zeros(nn,m-nn); zeros(m-nn,nn) zeros(m-nn)];
    E=U*Z*U';
end

W=eye(m)+E;

[M,p]=cholcov2((W+W')/2);

if isnan(p) || p>0
    error('cholcov2() fails!');
end



