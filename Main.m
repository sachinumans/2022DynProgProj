clearvars -except D; close all; clc; tic
plt = "0";
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

% if plt=="plot"
%     d = -1.5:0.03:1.5;
%     dNeg = d(d<=0);
%     dPos = d(d>=0);
%     figure(); hold on
%     plot([Zmean Zmean], [0 pdf(Z,Zmean)], "k:")
%     fill([fliplr(dNeg) dNeg],[zeros(1, size(dNeg,2)) pdf(Z, dNeg)], 'r', 'facealpha',.5)
%     fill([fliplr(dPos) dPos],[zeros(1, size(dPos,2)) pdf(Z, dPos)], 'g', 'facealpha',.5)
%     text(-.25, 0.5, num2str(cdf(Z, 0), 4))
%     text(.15, 0.5, num2str(1-cdf(Z, 0), 4))
%     title("Power consumption derivative")
%     xlabel("Z")
%     ylabel("P(Z)")
% end
% 
% toc

%% Simulation
z = [0];
q = [D.Load(1)];
c = [0];
m = [0];

alpha_top = 0.96;
alpha_bot = 0.87;
u_min = 4*10^3;
u_max = 5*10^3;
p_peak  =1.5;
c_max = 7*10^3;
delta= 1;



for t = 1:length(D.Load)-1
z(t+1) =  D.Load(t+1)-D.Load(t);

q(t+1) = q(t)+ z(t+1);

% c(t+1) = c(t)+ alpha_top*delta*u(t)*(u(t)>=0) ...
%     + alpha_bot*delta*u(t)*(u(t)<0);

m(t+1) = max(m(t),q(t+1));
end

m = reshape(m,1440,366);
m = max(m);
G_T = m(end);


%% Plots
for day = 1:1
    LoadDay = reshape(D.Load,1440,ceil(length(D.Load)/1440));
    [LoadDayMax,LoadDayMaxInd] = max(LoadDay(:,day));
    figure
    plot(LoadDay(:,day));
    hold on
    plot(LoadDayMaxInd,LoadDayMax,'ro');
    grid on
    xlabel("Time [min]");
    ylabel("Power[kW]");
    title("Day " + day);
    legend("Consumed Power","Maximum consumed power")
end