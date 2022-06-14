function [datCell, oneYear, ddt_oneYear, M, S] = ExtractGivenData(D, plt)

N = size(D,1);
idx_start = 1; % index of day start
dayCount = 1; % Day index
dayArr = convertCharsToStrings(D{:, 1});
dayStartIdx = cell([], 1); % Datacell with a cell for each day containing the starting index of that day
dayEndIdx = cell([], 1);

oneYear = table2array(D(:, 3));

datCell = cell([], 1); % Datacell with a cell for each day

if plt=="plot"; figure(); hold on; title("All days on top of each other"); end

while idx_start <= N
    day = convertCharsToStrings(D{idx_start,1});
    idx_end = find(dayArr == day, 1, 'last' ); % index of day end
    
    Y = table2array(D(idx_start:idx_end, 3));
    
    if plt=="plot"; plot(Y); end
    
    datCell{dayCount} = Y;
    dayStartIdx{dayCount} = idx_start;
    dayEndIdx{dayCount} = idx_end;
    
    idx_start = idx_end+1;
    
    dayCount = dayCount +1;
end
if plt=="plot"; hold off; end
% save datCell.mat datCell oneYear

%% Plot derivative of each day
n = max(size(datCell{1}));
X = zeros(max(size(datCell)), n);
ddt_oneYear = diff(oneYear);

if plt=="plot"
    figure(); hold on; 
    title("The derivative of each day on top of each other")
end
for d = 1:(dayCount-1)
    idxStart = dayStartIdx{d};
    idxEnd = dayEndIdx{d};
    if idxEnd == N
        X(d, 1:end-1) = ddt_oneYear(idxStart:N-1);
    else
        X(d, :) = ddt_oneYear(idxStart:idxEnd);
    end
    if plt=="plot"; plot(ddt_oneYear(idxStart:min(idxEnd,N-1))); end
    
end
if plt=="plot"; hold off; end

%% std of each minute of the day
S = std(X);

%% mean of each minute of the day
M = mean(X, 1);

%% Plots
if plt=="plot"
% Plot the entire year
figure(); hold on
title("The entire year")
days = (1:N).*365./N;
plot(days, oneYear);
hold off

% Plot derivative of the entire year
figure(); hold on
title("The derivate of entire year")
plot(days(2:end), ddt_oneYear);
hold off

% Plot mean of deriv
figure(); hold on
title("The mean of the derivative of each minute in the day")
plot(M);
hold off

% Plot std of deriv
figure(); hold on
title("The standard deviation of the derivative of each minute in the day")
plot(S);
hold off

end
end
