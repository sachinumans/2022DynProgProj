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


m(t+1) = max(m(t),q(t+1));
end

m = reshape(m,1440,366);
m = max(m);
G_T = m(end);


Z1 = cdf(Z, (-9.75:0.5:10.25)) - cdf(Z, (-10.25:0.5:9.75));

time = 60;
states = [20 20 14]; % current max (0.5kW steps); current consumption (0.5kW steps); current battery state(0.5kW steps)

 steps = 24; %hours

 U = zeros(states(1),states(2),states(3),steps);
 J = zeros(states(1),states(2),states(3),steps);

 max_discharge_rate = 5*10^3;
 max_charge_rate = 4*10^3;
 alpha_bot = 0.87;
 alpha_top = 0.96;

 term_cost_battery = 1;
 running_cost = 0;
 terminal_cost = 5;

 for term_cost_battery = 1:50
     for running_cost = 1:15
         for terminal_cost = 1:30
 
             for i = 1:states(1)
                for j = 1:states(2)
                    for k =1:states(3)
                      J(i,j,k,24) = terminal_cost*(i^2) + term_cost_battery*(14-k)^2;
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
                for j = 0:23       
            
                    u(i*24+j+1) = U(ceil(cost_max(1,i+1)), ceil(2*D.Load(i*24+j*60+1)), ceil(b(i*24+j+1)), j+1);
                    if j==23
                     u(i*24+j+1) = 9;   
                    end
                    if u(i*24+j+1) > 9           
                    b(i*24+j+2) =  max(b(i*24+j+1) - alpha_bot* (u(i*24+j+1)-9), 1);
                    end
            
                    if u(i*24+j+1) < 9 
                       b(i*24+j+2) = min(b(i*24+j+1) + alpha_top * (9- u(i*24+j+1)), 14);
                    end
            
                    if u(i*24+j+1) == 9 
                       b(i*24+j+2) = b(i*24+j+1);
                    end
            
                    for k=1:60
                            cost_max(1,i+1) = max(cost_max(1,i+1), D.Load(i*24+(j*60)+k) - (b(i*24+j+2)-b(i*24+j+1)/2));
                            cost_max(2,i+1) = max(cost_max(2,i+1), D.Load(i*24+(j*60)+k));
                    end
            
                end
            end
            
            price = [sum(cost_max(1,:))*1.5 sum(cost_max(2,:))*1.5];

            stuff(term_cost_battery,running_cost,terminal_cost) = price(1);
            
         end
     end
 end


%% Plots

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
for day = 1:1
    LoadDay = reshape(D.Load,1440,ceil(length(D.Load)/1440));
    [LoadDayMax,LoadDayMaxInd] = max(LoadDay(:,day));
    figure
    plot(LoadDay(:,day));
    hold on
    plot(LoadDayMaxInd,LoadDayMax,'ro');
    hold on
    plot(grid_power(:,day));
    hold on
    plot(max_ind(day),max_grid(day),'go');
    grid on
    xlabel("Time [min]");
    ylabel("Power[kW]");
    title("Day " + day);
    legend("Consumed Power","Maximum consumed power","Grid Power Consumed","Maximum grid power consumed")

    figure
    plot(b_power(:,day))
    hold on
    yline(umin,'r');
    yline(umax,'g');
    grid on
    xlabel("Time [min]");
    ylabel("Power[kW]");
    title("Day " + day);
    ylim([-7,7])
    legend("Battery change in power","Discharging limit","Charging limit")

    figure
    plot(b(:,day)./2)
    hold on
    yline(7);
    yline(0);
    grid on
    xlabel("Time [min]");
    ylabel("Capacity[kWh]");
    title("Day " + day);
    legend("Battery change in power","Discharging limit","Charging limit")
    ylim([-2,9])
end




