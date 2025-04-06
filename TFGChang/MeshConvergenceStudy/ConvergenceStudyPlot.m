%% Calculate E smooth and Absolute Relative Error (Case with FiltreLump) 

openfig('E_N_NACA9810_A0A10_PF_FL.fig');

h = findobj(gca, 'Type', 'line');
nx = get(h, 'XData'); 
E = get(h, 'YData');

windowSize = 5;  % Window size for the moving average
ySmooth = movmean(E, windowSize);

for i = 2:length(ySmooth)
    errorRelative(i-1) = abs((ySmooth(i) - ySmooth(i-1)) / ySmooth(i-1)) ;
end


%% Plot Original E & E_smooth vs nx and Absolute Relative Error vs nx (Case with FiltreLump) 

figure;

subplot(2,1,1);
plot(nx, E, 'b-', 'DisplayName', 'Original Data'); hold on;
plot(nx, ySmooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Moving Average');
legend; grid on;
title('Original and Smoothed Data');

subplot(2,1,2);
plot(nx(2:end), errorRelative, 'm-', 'LineWidth', 2, 'DisplayName', 'Absolute Relative Error');
legend; grid on;
title('Absolute Relative Error between Consecutive y_{smooth}');
ylim([0,0.1]);

%% Calculate E smooth and Absolute Relative Error (Case without FiltreLump) 

data = load("E_N9810_AoA10_PF_NFL.mat");

windowSize = 5;
ySmooth = movmean(data.E, windowSize);

for i = 2:length(ySmooth)
    errorRelative(i-1) = abs((ySmooth(i) - ySmooth(i-1)) / ySmooth(i-1)) ;
end

%% Plot Original E & E_smooth vs nx and Absolute Relative Error vs nx (Case without FiltreLump) 

figure;

subplot(2,1,1);
plot(nx, data.E, 'b-', 'DisplayName', 'Original Data'); hold on;
plot(nx, ySmooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Moving Average');
ylim([0,1.4]);
legend; grid on;
title('Original and Smoothed Data');

subplot(2,1,2);
plot(nx(2:end), errorRelative, 'm-', 'LineWidth', 2, 'DisplayName', 'Absolute Relative Error');
legend; grid on;
title('Absolute Relative Error between Consecutive y_{smooth}');
ylim([0,0.1]);