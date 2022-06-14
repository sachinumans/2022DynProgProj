clearvars -except D; close all; clc; tic
plt = "plot";
%% Define vars
cmax = 7; %kWh
ubAlpha = 0.96;
lbAlpha = 0.87;
umin = -4; %kWh
umax = 5; %kWh
ppeak = 1.50; % euro

s=whos;
vars = cell2struct({s.name}.',{s.name});

%% Load cvs
if ~exist("D", "var"); D = readtable("data.csv"); end

[datCell, oneYear, ddt_oneYear, M, S] = ExtractGivenData(D, plt);

%% Determine characteristics of randomness
Zmean = mean(ddt_oneYear);
Zstd = std(ddt_oneYear);

Z = makedist('Normal','mu',Zmean,'sigma',Zstd);

if plt=="plot"
    d = -1.5:0.03:1.5;
    dNeg = d(d<=0);
    dPos = d(d>=0);
    figure(); hold on
    plot([Zmean Zmean], [0 pdf(Z,Zmean)], "k:")
    fill([fliplr(dNeg) dNeg],[zeros(1, size(dNeg,2)) pdf(Z, dNeg)], 'r', 'facealpha',.5)
    fill([fliplr(dPos) dPos],[zeros(1, size(dPos,2)) pdf(Z, dPos)], 'g', 'facealpha',.5)
    text(-.25, 0.5, num2str(cdf(Z, 0), 4))
    text(.15, 0.5, num2str(1-cdf(Z, 0), 4))
    title("Power consumption derivative")
    xlabel("Z")
    ylabel("P(Z)")
end

toc