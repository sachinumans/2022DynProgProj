% clear all; 
close all; clc;
% Load cvs
if ~exist("D", "var"); D = readtable("data.csv"); end
N = size(D,1);
idx_start = 1; % index of day start
dayCount = 1; % Day index
dayArr = convertCharsToStrings(D{:, 1});

oneYear = table2array(D(:, 3));

datCell = cell([], 1); % Datacell with a cell for each day

figure(); hold on;

while idx_start <= N
    day = convertCharsToStrings(D{idx_start,1});
    idx_end = find(dayArr == day, 1, 'last' ); % index of day end
    
    Y = table2array(D(idx_start:idx_end, 3));
    plot(Y)
    
    datCell{dayCount} = Y;
    
    idx_start = idx_end+1;
    
    dayCount = dayCount +1;
end

save datCell.mat datCell oneYear