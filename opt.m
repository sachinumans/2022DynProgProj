clear all
close all
clc

if ~exist("D", "var"); D = readtable("data.csv"); end

load('Z.mat');

time = 60;
states = [20 20 14];

steps = 24;

 U = zeros(states(1),states(2),states(3),steps);
 J = zeros(states(1),states(2),states(3),steps);

 max_consumption = 5;
 term_cost_battery = 1;
 running_cost = 0;
 
 for i = 1:states(1)
    for j = 1:states(2)
        for k =1:states(3)
          J(i,j,k,1) = max_consumption*(i^2) + term_cost_battery*(14-k)^2;
        end
    end 
 end

for p = 2: steps
     for i = 1:states(1)
        for j = 1:states(2)
            for k =1:states(3)
              m = ones(17,1)*(10^9);
                
                for n = 8:1
                    if j <(states(2) - n) && k < (states(3)-n)
%                         m(9-n) = random(Z,[1 20])*J(:,j,k+n+1,p-1);
                          Z1(:) = squeeze(Z(i,j-n, :));
                          m(9-n) = Z1 * J(:,j,k+n,p-1);
                    end
                end
%                 m(9) = random(Z,[1 20])*J(:,j,k,p-1);
                Z1(:) = squeeze(Z(i,j, :));
                m(9) = Z1*J(:,j,k,p-1);
                for n = 1:8
                    if (j > (n-1)) && (k > (n-1))
%                       m(9+n) = random(Z,[1 20])*J(:,j,k-n+1,p-1);  
                        Z1(:) = squeeze(Z(i,j, :));
                        m(9+n) = Z1*J(:,j,k,p-1);
                    end               
                end
                [J(i,j,k,p),U(i,j,k,p)] = min(m);                
            end
        end 
     end
end


%% Implementation
alpha_bot = 0.87;
alpha_top = 0.96;

cost_max = ones(2, 366);
c = ones(1,8785);
b = ones(1,8785);
u = ones(1,8785);

for i=0:365
    for j = 0:23
        c(i*24+j+1) =  D.Load(i*24+j*60+1); 
    
        u(i*24+j+2) = U(cost_max(1,i+1), floor(2*c(i*24+j+1)), floor(b(i*24+j+1)), 24-j);

        if u(i*24+j+2) > 9            
        u(i*24+j+2) =  b(i*24+j+1) - alpha_bot* (u(i*24+j+2)-9);
        end
        if u(j+2) < 9 
           b(i*24+j+2) = b(i*24+j+1) - alpha_top * (u(i*24+j+2)-9);
        end

        for k=1:60
                cost_max(1,i+1) = max(cost_max(1,i+1), D.Load(i*24+(j*60)+k) - u(j+2));
                cost_max(2,i+1) = max(cost_max(2,i+1), D.Load(i*24+(j*60)+k));
        end

    end
end

price = [sum(cost_max(1,:))*1.5 sum(cost_max(2,:))*1.5];
