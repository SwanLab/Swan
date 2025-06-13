%% Superposition of Linear Hardening and Softening Laws
clc;
clear;
close all;

% Material parameters
E0 = 210;
r0 = 10;
r1 = 20;
A = 1;


% Applied displacement values
uVals = 0:1e-3:0.2;

% Define hardening and softening parameters
q_infVals = [10, -0.1]; % H = +0.5 (hardening), H = -0.5 (softening)
H_vals = [0.5, -0.5];
labels = {'Hardening (q_\infty = 0.5)', 'Softening (q_\infty = -0.5)'};

% Initialize storage for both cases
dmg_all = {};
q_all = {};
F_all = {};
tau_all = {};

for h = 1:length(q_infVals)
    qInf = q_infVals(h);
    H = H_vals(h)
    rOld = r0;

    dmg = [];
    qVec = [];
    F = [];
    tauVec = [];

    for i = 1:length(uVals)
        u = uVals(i);

        % Compute tau (stress measure)
        tau = sqrt(u * E0 * u);

        % Update internal variable r
        if tau > rOld
            r = tau;
        else
            r = rOld;
        end

        % Compute hardening/softening function q
        if r >= r1
            q = qInf;
        else
            q =  qInf - (qInf - r0)*exp(A*(1-r/r0));
        end

        % Damage calculation
        d = 1 - q / r;
        d = min(max(d, 0), 1);  % Clamp d to [0,1]

        % Compute reaction force
        reac = (1 - d) * E0 * u;

        % Update rOld
        rOld = r;

        % Store results
        dmg(end+1) = d;
        qVec(end+1) = q;
        F(end+1) = reac;
        tauVec(end+1) = tau;
    end

    % Store for comparison
    dmg_all{end+1} = dmg;
    q_all{end+1} = qVec;
    F_all{end+1} = F;
    tau_all{end+1} = tauVec;
end

% Plot Damage vs Displacement
figure();
hold on;
for i = 1:2
    plot(tau_all{i}, dmg_all{i}, 'LineWidth', 2);
end
xlabel('Damage evolution parameter (r)', 'Interpreter', 'latex');
ylabel('Damage (d)', 'Interpreter', 'latex');
title('Damage variable vs damage evolution parameter');
legend(labels, 'Location', 'best');
grid on;

% Plot q vs Tau
figure();
hold on;
for i = 1:2
    plot(tau_all{i}, q_all{i}, 'LineWidth', 2);
end
xlabel('Damage evolution parameter (r)', 'Interpreter', 'latex');
ylabel('Hardening law (q)', 'Interpreter', 'latex');
title('Exponential hardening law as a function of the damage evolution parameter');
legend(labels, 'Location', 'best');
grid on;

% Plot Force vs Displacement
figure();
hold on;
for i = 1:2
    plot(tau_all{i}, F_all{i}, 'LineWidth', 2);
end
xlabel('Damage evolution parameter (r)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
title('Reaction force as a function of the damage evolution parameter');
legend(labels, 'Location', 'best');
grid on;
% Plot Displacement vs r
figure();
hold on;
for i = 1:2
    plot(uVals,max(tau_all{i},r0), 'LineWidth', 2);
end
xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Damage evolution parameter (r)', 'Interpreter', 'latex');
title('Damage evolution parameter as a function of displacement');
legend(labels, 'Location', 'best');
grid on;