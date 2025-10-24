clear; close all; clc;


%% 1. Parameters

desired_std = 1;
N = 64e3;
fs = 8e3;

f0 = fs/N;

disp('Experiment duration (s):');
disp(N/fs);

%% 2. Signal generation

u0 = desired_std * randn(N,1); 

if ~exist('signals','dir')
    mkdir('signals');
end
save(fullfile('signals','task_1_1.mat'),'u0');

