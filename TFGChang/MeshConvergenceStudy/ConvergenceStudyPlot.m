%% Calculate E smooth and Absolute Relative Error (Case with FiltreLump, Stokes)

openfig('E_N_NACA9810_A0A10_PF_FL.fig');

h = findobj(gca, 'Type', 'line');
nx = get(h, 'XData');
E = get(h, 'YData');

windowSize = 5;  % Window size for the moving average
ySmooth = movmean(E, windowSize);

for i = 2:length(ySmooth)
    errorRelative(i-1) = abs((ySmooth(i) - ySmooth(i-1)) / ySmooth(i-1)) ;
end


%% Plot Original E & E_smooth vs nx and Absolute Relative Error vs nx (Case with FiltreLump, Stokes)

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

%% Calculate E smooth and Absolute Relative Error (Case without FiltreLump, Stokes)

data = load("E_N9810_AoA10_PF_NFL.mat");

windowSize = 5;
ySmooth = movmean(data.E, windowSize);

for i = 2:length(ySmooth)
    errorRelative(i-1) = abs((ySmooth(i) - ySmooth(i-1)) / ySmooth(i-1)) ;
end

%% Plot Original E & E_smooth vs nx and Absolute Relative Error vs nx (Case without FiltreLump, Stokes)

figure;

%subplot(2,1,1);
plot(nx, data.E, 'b-', 'DisplayName', 'Original Data'); hold on;
plot(nx, ySmooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Moving Average');
%ylim([0,1.4]);
legend; grid on;
title('Original and Smoothed Data');
%
subplot(2,1,2);
plot(nx(2:end), errorRelative, 'm-', 'LineWidth', 2, 'DisplayName', 'Absolute Relative Error');
legend; grid on;
title('Absolute Relative Error between Consecutive y_{smooth}');
ylim([0,0.1]);


%% Plot E and Absolute Relative Increment vs nx (Navier-Stokes)

data = readmatrix("E_nx_NS.txt");
nx = 20:20: size(data,1)*20;

E = data(:,end);

for i = 2:length(E)
    errorRelative(i-1) = abs((E(i) - E(i-1)) / E(i-1)) ;
end

figure;

plot(nx,data(:,end));
grid on;
title('E vs. n_x');

figure;
plot(nx(2:end), errorRelative);
grid on;
title('Absolute Relative Error between consecutive E vs. n_x');
ylim([0,0.15]);