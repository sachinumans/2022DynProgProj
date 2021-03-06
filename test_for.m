clearvars -except D; close all; clc; tic
plt = "plodt";
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


Z1 = cdf(Z, (1:20)./20.*2-1) - cdf(Z, (0:19)./20.*2-1);


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

%% Simulation
z = [0];
q = [D.Load(1)];
c = [0];
m = [0];
cost_peak  =1.5;
c_max = 7;
delta= 1;



for t = 1:length(D.Load)-1
z(t+1) =  D.Load(t+1)-D.Load(t);

q(t+1) = q(t)+ z(t+1);


m(t+1) = max(m(t),q(t+1));
end

m = reshape(m,1440,366);
m = max(m);
G_T = m(end);


Z1 = cdf(Z, (-9.75:0.5:10.25)) - cdf(Z, (-10.25:0.5:9.75));

time = 60;
states = [20 20 14]; % current max (0.5kW steps); current consumption (0.5kW steps); current battery state(0.5kW steps)

 steps = 24*time; %minutes

 U = zeros(states(1),states(2),states(3),steps);
 J = zeros(states(1),states(2),states(3),steps);

 alpha_bot = 0.87;
 alpha_top = 0.96;

 term_cost_battery = 1000;
 running_cost = 0;
 terminal_cost = 500;

 
 for i = 1:states(1)
    for j = 1:states(2)
        for k =1:states(3)
          J(i,j,k,steps) = terminal_cost*(i^2) + term_cost_battery*(14-k)^2;
        end
    end 
 end

for p = steps-1:-1:1
     for i = 1:states(1)
        for j = 1:states(2)
            for k =1:states(3)
              m = ones(19,1)*(10^9);
              Zloc = Z1(22-j:21-j+20);
              Zloc(1) = Zloc(1) + sum(Z1(1:22-j-1));
              Zloc(end) = Zloc(end) + sum(Z1(21-j+20+1:end));
                for n = 1:8                   
                   if (j+n <= 20) && (k+n <= 14)
                        m(n) = running_cost*(14 - k)*2 + Zloc*J(:,j+n,k+n,p+1); 
                   end
                end              
              
                m(9) = running_cost*(14 - k)*2 + Zloc*J(:,j,k,p+1);
                
                for n = 1:10
                    if (j-n >= 1) && (k-n >=1)
                        m(9+n) = running_cost*(14 - k)*2 + Zloc*J(:,j-n,k-n,p+1); 
                    end               
                end

                [J(i,j,k,p),U(i,j,k,p)] = min(m);

            end
        end 
     end
end

%% Implementation

cost_max = ones(2, 366);
b = ones(1,8785);
u = ones(1,8785)*9;

for i=0:365
    for j = 0:steps-1       

        u(i*steps+j+1) = U(floor(cost_max(1,i+1)), ceil(2*min(D.Load(i*steps+j+1),10)), floor(b(i*steps+j+1)), j+1);
        if j==steps-1
         u(i*steps+j+1) = 9;   
        end
        if u(i*steps+j+1) > 9           
        b(i*steps+j+2) =  max(b(i*steps+j+1) - alpha_bot* (u(i*steps+j+1)-9), 1);
        end

        if u(i*steps+j+1) < 9 
           b(i*steps+j+2) = min(b(i*steps+j+1) + alpha_top * (9- u(i*steps+j+1)), 14);
        end

        if u(i*steps+j+1) == 9 
           b(i*steps+j+2) = b(i*steps+j+1);
        end


        cost_max(1,i+1) = max(cost_max(1,i+1), D.Load(i*steps+j+1) - (b(i*steps+j+2)-b(i*steps+j+1))/2);
        cost_max(2,i+1) = max(cost_max(2,i+1), D.Load(i*steps+j+1));

    end
end

price = [sum(cost_max(1,:))*1.5 sum(cost_max(2,:))*1.5];


%% Plots

% Battery Power Consumption
b_power = diff(b)./2;
% Grid Power
for i=0:365
    for j = 0:steps-1  
        grid_power(i*steps+j+1) = max(D.Load(i*steps+j+1)+b_power(i*steps+j+1),0);        
    end
end

grid_power = reshape(grid_power,steps,366);
[max_grid,max_ind] = max(grid_power);

b_power = reshape(b_power,steps,366);

b = reshape(b(1:end-1),steps,366);
%% Plots
for day = 1:30
    LoadDay = reshape(D.Load,1440,ceil(length(D.Load)/1440));
    bat = LoadDay(:,day) - grid_power(:,day);
    [LoadDayMax,LoadDayMaxInd] = max(LoadDay(:,day));
    figure
    subplot(3,1,1)
    plot(LoadDay(:,day),'r','LineWidth',1.5);
    hold on
    plot(LoadDayMaxInd,LoadDayMax,'ro','LineWidth',1.5);
    hold on
    plot(grid_power(:,day),'b','LineWidth',1.5);
    hold on
    plot(max_ind(day),max_grid(day),'bo','LineWidth',1.5);
    grid on
    xlabel("Time [min]");
    ylabel("Power[kW]");
    title("Day " + day);
    legend("Consumed Power","Maximum consumed power","Grid Power Consumed","Maximum grid power consumed")

    subplot(3,1,2)
    plot(bat,'r','LineWidth',1.5)
    hold on
    yline(umin,'b','LineWidth',1.5);
    yline(umax,'k','LineWidth',1.5);
    grid on
    xlabel("Time [min]");
    ylabel("Power[kW]");
    title("Day " + day);
    ylim([-7,7])
    legend("Battery change in power","Discharging limit","Charging limit")

    subplot(3,1,3)
    plot(b(:,day)./2,'r','LineWidth',1.5)
    hold on
    yline(7,'k','LineWidth',1.5);
    yline(0,'b','LineWidth',1.5);
    grid on
    xlabel("Time [min]");
    ylabel("Capacity[kWh]");
    title("Day " + day);
    legend("Battery change in power","Discharging limit","Charging limit")
    ylim([-2,9])
end




