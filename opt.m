clear all
close all
clc

if ~exist("D", "var"); D = readtable("data.csv"); end

load('matlab.mat');
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
%                 cost_max(1,i+1) = max(cost_max(1,i+1), D.Load(i*24+(j*60)+k) - 0.5*(u(i*24+j+1)-9));
                cost_max(1,i+1) = max(cost_max(1,i+1), D.Load(i*24+(j*60)+k) - (b(i*24+j+2)-b(i*24+j+1)/2));
                cost_max(2,i+1) = max(cost_max(2,i+1), D.Load(i*24+(j*60)+k));
        end

    end
end

price = [sum(cost_max(1,:))*1.5 sum(cost_max(2,:))*1.5];


%% Plots

% Battery Power Consumption
% b_power = kron(diff(b),ones(1, 60))./2;
b_power = diff(b)./2;
% Grid Power
for i=0:365
    for j = 0:23  
        grid_power(i*24+j+1) = max(D.Load(i*24+j*60+1)-b_power(i*24+j+1),0);        
    end
end

max_grid = reshape(grid_power,24,366);
max_grid = max(max_grid);








