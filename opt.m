clear all
close all
clc

if ~exist("D", "var"); D = readtable("data.csv"); end

load('matlab.mat');

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
                        m(9-n) = random(Z,[1 20])*J(:,j,k+n+1,p-1);
                    end
                end
                m(9) = random(Z,[1 20])*J(:,j,k,p-1);
                for n = 1:8
                    if (j > (n-1)) && (k > (n-1))
                      m(9+n) = random(Z,[1 20])*J(:,j,k-n+1,p-1);   
                    end               
                end
                [J(i,j,k,p),U(i,j,k,p)] = min(m);                
            end
        end 
     end
end