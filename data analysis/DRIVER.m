%% Predicting the Frequencies of drug side effects
%  Diego Galeano, Shantao Li, Mark Gerstein, Alberto Paccanaro
%  Nature Communications
%
%  Instructions
%  ------------
%
%  This file contains code to run the decomposition algorithm and generate 
%  predictions for a given drug-side effect of interes.
%
%  Code Copyright (C) 2020 by Diego Galeano.

%% Initialization
% Clears the variables in the environment and add folders to the path
clc; clear all; close all;
addpath('data\'); addpath('source\');


%% ================== Section 1: Load Dataset  =======================
load('FrequencyData.mat');


%% ================== Section 2: Train the model  ====================
% Our matrix decomposition model parameters
k = 10;       % number of signature components 
alpha = 0.05; % confidence on the zeros

tic;
fprintf('Training the model...\n\n');
[ W, H ] = DecompositionAlgorithm( R, k, alpha );
toc;

% prediction model
Rhat = W*H; % drug signatures x side effect signatures

%% ================== Section 3: Predicting the specific frequency class  ===============
classes = {'zero', 'very rare', 'rare', 'infrequent', 'frequent', 'very frequent'};

mydrug = 'alfentanil';
mysideeffect = 'nervousness';

idx =  strcmp(drugs, mydrug);
idy = strcmp(sideeffects, mysideeffect);

fprintf('-------------drug: %s, side effect: %s --------------\n\n', mydrug, mysideeffect)
fprintf('Original entry value for the pair %s\n', upper(classes{R(idx, idy) + 1}));
fprintf('Predicted frequency class for the pair %s\n', upper(classes{predictedFrequencyClass(Rhat(idx, idy))+1}));

%
%
%
%  Code Copyright (C) 2020 by Diego Galeano.
%  Any doubt write to the corresponding author alberto.paccanaro@rhul.ac.uk
%  or to the first author diego.galeano@fgv.br (https://diegogalpy.github.io/)
