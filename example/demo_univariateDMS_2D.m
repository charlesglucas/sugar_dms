
%This code shows how to use the toolbox SUGAR_DMS minimizing the Discrete
%Mummford-Shah functional:
%           1/2*||z-u||^2 + beta*||(1-e).*Du||^2 + lambda*R(e)
%with automatic tuning of hyperparamters
%
% Implementation C.G. LUCAS, ENS Lyon
clc
clear all
close all

format compact
addpath(genpath('./'))

%% noisy RGB image
im = imread('data/ellipse-128.pgm');
x = double(im)/255;
[n1,n2] = size(x);
stdn = 0.01; % noise level
z = x + stdn*randn(n1,n2); z = double(z);

%% DMS with automatic selection of parameters using BFGS
param = struct;
param.R = 5; % number of realizations of the Monte Carlo vector
param.sigma = stdn; % noise level estimation if empty

%% SUGAR DMS
[Lambda,crit] = bfgs_sugar_dms(z,param);
[u,e,~] = DMS_2D(z,Lambda(1),Lambda(2));

%% Plots DMS estimates

figure(1); clf;
subplot(1,3,1);
imshow(x); axis on; title('Original image')
set(gca,'xtick',[]); set(gca,'ytick',[]);
subplot(1,3,2);
imshow(z); axis on; title('Noisy image')
set(gca,'xtick',[]); set(gca,'ytick',[]);
pbaspect([1 1 1])
subplot(1,3,3);
imshow(u); axis on; title('D-MS estimate')
set(gca,'xtick',[]); set(gca,'ytick',[]); hold on;
plot_contours(e,{'LineStyle','-','LineWidth',1,'Color','r'})
