clearvars -except D; close all; clc; tic
plt = "plodt";
%% Define vars
cmax = 7; %kWh
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


toc

%% Simulation
Z1 = cdf(Z, (-9.75:0.5:10.25)) - cdf(Z, (-10.25:0.5:9.75)); %draw from the identified distribution

states = [20 20 14]; % current max (0.5kW steps); current consumption (0.5kW steps); current battery state(0.5kW steps)

 steps = 24; %hours

 U = zeros(states(1),states(2),states(3),steps);
 J = zeros(states(1),states(2),states(3),steps);


 alpha_bot = 0.87;
 alpha_top = 0.96;

 term_cost_battery = 5;
 running_cost = 0;
 terminal_cost = 5;

 
 for i = 1:states(1)
    for j = 1:states(2)
        for k =1:states(3)
          J(i,j,k,24) = terminal_cost*(i^2) + term_cost_battery*(14-k)^2;
        end
    end 
 end

% battery_charge_state(1) = battery will be fully charged
% battery_charge_state(1-8) battery will be charged in steps of 0.5kWh till the charging rate
% limit of 4kWh is reached in battery_charge_state(8)
% battery_charge_state(9) - battery is neither charged nor discharged
% battery_charge_state(10-20) battery will be discharged in steps of 0.5kWh till the discharging rate
% limit of 5kWh is reached in battery_charge_state(20)

for p = steps-1:-1:1 % iterate through each hour
     for i = 1:states(1) %iterate through each current max
        for j = 1:states(2) %iterate through each current consumption
            for k =1:states(3) %iterate through each battery level
              battery_charge_state = ones(19,1)*(10^9);
              Zloc = Z1(22-j:21-j+20);
              Zloc(1) = Zloc(1) + sum(Z1(1:22-j-1));
              Zloc(end) = Zloc(end) + sum(Z1(21-j+20+1:end));
                for n = 1:8                   
                   if (j+n <= 20) && (k+n <= 14)
                        battery_charge_state(n) = running_cost*(14 - k)*2 + Zloc*J(:,j+n,k+n,p+1); % candidate cost function, the stochasticity Zloc is added to J of t+1 
                   end
                end              
              
                battery_charge_state(9) = running_cost*(14 - k)*2 + Zloc*J(:,j,k,p+1); % candidate cost function, the stochasticity Zloc is added to J of t+1
                
                for n = 1:10
                    if (j-n >= 1) && (k-n >=1)
                        battery_charge_state(9+n) = running_cost*(14 - k)*2 + Zloc*J(:,j-n,k-n,p+1); % candidate cost function, the stochasticity Zloc is added to J of t+1
                    end               
                end

                [J(i,j,k,p),U(i,j,k,p)] = min(battery_charge_state); % find the U that minimizes the cost function J

            end
        end 
     end
end

%% Implementation

cost_max = ones(2, 366);
b = ones(1,8785);
u = ones(1,8785)*9;

for i=0:365 %iterate through the whole year
    for j = 0:23  %iterate through every hour     

        u(i*24+j+1) = U(ceil(cost_max(1,i+1)), ceil(2*D.Load(i*24+j*60+1)), ceil(b(i*24+j+1)), j+1); %choose optimal input based on the current state and data
        if j==23
         u(i*24+j+1) = 9;   
        end
        if u(i*24+j+1) > 9           
        b(i*24+j+2) =  max(b(i*24+j+1) - alpha_bot* (u(i*24+j+1)-9), 1); %discharge battery according to optimal input
        end

        if u(i*24+j+1) < 9 
           b(i*24+j+2) = min(b(i*24+j+1) + alpha_top * (9- u(i*24+j+1)), 14); %charge battery according to optimal input
        end

        if u(i*24+j+1) == 9 
           b(i*24+j+2) = b(i*24+j+1); %keep battery charge according to optimal input
        end
        %Calculate maximum cost per respective day
        for k=1:60
                cost_max(1,i+1) = max(cost_max(1,i+1), (D.Load(i*24+(j*60)+k) - (b(i*24+j+2)/2-b(i*24+j+1)/2))/2);
                cost_max(2,i+1) = max(cost_max(2,i+1), D.Load(i*24+(j*60)+k));
        end

    end
end
% Price for the whole year
price = [sum(cost_max(1,:))*1.5 sum(cost_max(2,:))*1.5];
disp("The total price for the whole year without control is:")
price(2)
disp("The total price for the whole year withcontrol is:")
price(1)
disp("The price was imporoved with")
percent = 100-(price(1)/price(2))*100;
disp(percent + "%")

% Battery Power Consumption
b_power = diff(b)./2;
% Grid Power
for i=0:365
    for j = 0:23  
        grid_power(i*24+j+1) = max(D.Load(i*24+j*60+1)-b_power(i*24+j+1),0);        
    end
end

grid_power = reshape(grid_power,24,366);
grid_power = kron(grid_power,ones(60, 1));
[max_grid,max_ind] = max(grid_power);

b_power = reshape(b_power,24,366);
b_power = kron(b_power,ones(60, 1));

b = reshape(b(1:end-1),24,366);
b = kron(b,ones(60, 1));



%% Plots
for day = 1:2
    LoadDay = reshape(D.Load,1440,ceil(length(D.Load)/1440));
    [LoadDayMax,LoadDayMaxInd] = max(LoadDay(:,day));
    figure
    subplot(2,1,1)
    bat = -(LoadDay(:,day) - grid_power(:,day));
    bat(bat<-4) = -4;
    bat(bat>5) = 5;
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

    subplot(2,1,2)
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

%     subplot(3,1,3)
%     plot(b(:,day)./2,'r','LineWidth',1.5)
%     hold on
%     yline(7,'k','LineWidth',1.5);
%     yline(0,'b','LineWidth',1.5);
%     grid on
%     xlabel("Time [min]");
%     ylabel("Capacity[kWh]");
%     title("Day " + day);
%     legend("Battery change in power","Discharging limit","Charging limit")
%     ylim([-2,9])
end




