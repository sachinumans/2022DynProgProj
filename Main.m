clear all; close all; clc;
%% Define vars

cmax = 7; %kWh
ubAlpha = 0.96;
lbAlpha = 0.87;
umin = -4; %kWh
umax = 5; %kWh
ppeak = 1.50; % euro

s=whos;
vars = cell2struct({s.name}.',{s.name});

%% Load data
load datCell

%%
