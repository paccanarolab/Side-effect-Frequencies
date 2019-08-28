%% Predicting the Frequencies of drug side effects
%  Diego Galeano, Shantao Li, Mark Gerstein, Alberto Paccanaro
%
%  Instructions
%  ------------
%
%  This file contains code that reproduce the prediction results on the 
%  held-out and post-marketing test sets. 
%
%  Copyright (C) 2019 by Diego Galeano.

%% Initialization
% Clears the variables in the environment and add folders to the path
clc; clear all; close all;
addpath('data\'); addpath('source\');


%% ================== Section 1: Load Dataset  ===================
%  We start by loading three matrices: training (R_train), and two test sets used
%  in the paper: R_TestHoldOut and R_TestPosMarket. 
%  It also loads two cell arrays containing the drug generic 
%  names (drugNames) and side effect MedDRA Preferred Terms (SELabels) corresponding to the 
%  rows and columns of the matrices, respectively.
%  Load the datasets
load('data.mat');


%% ================== Section 2: Train the model  ====================
% Our matrix decomposition model parameters
K = 10;       % number of latent features
alpha = 0.05; % confidence on the zeros

tic;
fprintf('Training the model...\n\n');
[ W, H ] = DecompositionAlgorithm( R_Train, K, alpha );
toc;

% prediction model
Res = W*H; % drug signatures x side effect signatures

%% ================== Section 3: Evaluation in HeldOut  ===============
[AUC_heldOut] = getAUROC(Res, R_TestHoldOut, R_Train);
[RMSE_heldOut] = getRMSE(Res, R_TestHoldOut);

fprintf(' AUC %.3f and RMSE %.3f\n', AUC_heldOut, RMSE_heldOut);

%% ================== Section 4: Evaluation in PostMarket  ===============
[AUC_postMarket] = getAUROC(Res, R_TestPostMarket, R_Train);
[RMSE_postMarket] = getRMSE(Res, R_TestPostMarket);

fprintf(' AUC %.3f and RMSE %.3f\n', AUC_postMarket, RMSE_postMarket);

%% ================== Section 5: Reproduce Fig. 2  ====================
% Histogram of predicted scores for each of the classes in the HeldOut test
% set and the Post Market associations.

frequencyClasses = {'very rare', 'rare', 'infrequent', 'frequent', 'very frequent'};

% colors.
cmap = [  0,109,44;...
          254,217,118;...
          254,178,76;...
          253,141,60;...
          240,59,32;...
          189,0,38]./255;
      
figure(2);
subplot(6,1, 1);
histogram(Res(R_TestPostMarket > 0), 50,...
                'Normalization', 'pdf',...
                'FaceColor', cmap(1, :));
legend('Post Marketing');
ylabel('pdf');
xlabel('scores');
grid on;
xlim([0 max(Res(:))]);
ylim([0 1]);
    
for class = 1:5
    subplot(6,1, class + 1);
    histogram(Res(R_TestHoldOut == class),50,...
                                       'Normalization', 'pdf',...
                                       'FaceColor', cmap(class + 1, :));
    legend(frequencyClasses(class));
    ylabel('pdf');
    xlabel('scores');
    grid on;
    xlim([0 max(Res(:))]);
    ylim([0 1]);
end



