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


%%  This is a simple demo showing how to use the LCFV optimization
%   This script assumes you already have computed Fisher Vectors.
%   If you need to compute Fisher Vectors from images, please refer to:
%       run_synthetic_end2end_demo.m
%       run_cifar10_end2end_demo.m (uses a subset of CIFAR-10)
%   
%   In this demo, we simply show how to use the two core LCFV functions:
%       solve_LCFV1() and solve_LCFV2()

clear; clc; close all;

%% load data
if ~exist('../data/example_data.mat', 'file')
    error('Example data not found!');
end
load('../data/example_data.mat');

%% get label comparison matrix
C1 = repmat(labels, 1, length(labels));
C2 = repmat(labels', length(labels), 1);
C = double(C1 == C2);

%% compute LCFV
alpha = 10;
G = fv';

% LCFV1
tic;
[M1, W1] = solve_LCFV1(G, C, alpha);
t1 = toc;
LCFV1 = M1 * G;
fprintf('LCFV1 took %f seconds \n', t1);

% LCFV2
tic;
M2 = solve_LCFV2(G, C, alpha);
t2 = toc;
LCFV2 = M2 * G;
fprintf('LCFV2 took %f seconds \n', t2);

