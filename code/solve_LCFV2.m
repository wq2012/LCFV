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


%%  This is the implementation of LCFV2
%   Method: 
%       G'*M'*M*G=G'*G+C, and M=[I;E], and C=A'*A, 
%       then solve M*G=[G;A],
%       or solve E in E*G=A
%   G: m*n Fisher vectors
%   C: n*n label comparison matrix
%   alpha: parameter
%   M: m*m resulting transformation matrix to be applied on Fisher vectors
%   M*G is called LCFV
%   A: l*n matrix
%   E: l*m matrix
%   m: dimension of features
%   n: number of instances
%   l: the new dimensions introduced by A

function M=solve_LCFV2(G,C,alpha)

[m,n]=size(G);

C=C*alpha;

sv_t=0.01; % threshold for small singular values

A=cholcov(C);

l=size(A,1);

if n>=m % over-determined
    disp('LCFV2: over-determined problem');
    E = A / G;
else % under-determined
    disp('LCFV2: under-determined problem');
    [U,S,V] = svd(G); % X = U*S*V'
    nn=sum(diag(S)>sv_t); % remove small singular values
    S=S(:,1:nn); % remove small singular values
    V=V(:,1:nn); % remove small singular values
    A=A*V;
    Sinv=diag(1./diag(S));
    Z1=A*Sinv;
    Z=[Z1 zeros(l,m-nn)];
    E=Z*U';
end

M=[eye(m);E];


