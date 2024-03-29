% Clear the environment
clc; clear; close all;

% Add all subfolders to the path
addpath(genpath(pwd));

% Load the raw butterfly network
load('Raw_butterfly_network.mat');

% Run NR and NE on the raw network
W_butterfly_NR = NR(W_butterfly0, 0.5); % Apply Network Regularization
W_butterfly_NE_NR = Network_Enhancement(W_butterfly_NR); % Apply Network Enhancement

% Calculate and print the accuracy on the raw and processed networks
[~, acc_raw] = CalACC(W_butterfly0, labels); % Accuracy on the raw network
[~, acc_NE_NR] = CalACC(W_butterfly_NE_NR, labels); % Accuracy on the enhanced and regularized network

% Print the results
fprintf('The accuracy on raw network is %6.4f \n', acc_raw);
fprintf('The accuracy after applying NE+NR is %6.4f \n', acc_NE_NR);

