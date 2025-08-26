%% Superposition of Linear Hardening and Softening Laws
clc;
clear;
close all;clear;

% Material parameters
E0 = 210;
r0 = 10;
r1 = 2000000;
A = 1;


% Applied displacement values
uVals = 0:1e-3:10;

% Define hardening and softening parameters
q_infVals = [20, -0.1]; % H = +0.5 (hardening), H = -0.5 (softening)
H_vals = [0.5, -0.5];
labels = {'Hardening: q_{\infty} = 20', 'Softening: q_{\infty} = -0.1'};

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
            q = r0 + H*(r1-r0);
        else
            q =  qInf - (qInf - r0)*exp(A*(1-r/r0));
            %q = r0 + H*(r-r0);
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

load('1ELEM_TRACTION_LINEAR_HARDENING.mat','data');
dataHard = data;
load('1ELEM_TRACTION_LINEAR_SOFT.mat','data');
dataSoft = data;

% Plot Damage vs Displacement
figure();
hold on;

plot(uVals, dmg_all{1}, 'LineWidth', 2);
plot(uVals, dmg_all{2}, 'LineWidth', 2);

% plot(dataHard.displacement.value,dataHard.damage.maxValue,'o','MarkerSize',10,'MarkerIndices',1:30:length(dataHard.damage.maxValue),Color='#0072BD')
% plot(dataSoft.displacement.value,dataSoft.damage.maxValue,'square','MarkerSize',10,'MarkerIndices',1:20:length(dataSoft.damage.maxValue),Color='#D95319')

xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Damage (d)', 'Interpreter', 'latex');
% title('Damage variable vs damage evolution parameter');
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')

% Plot q vs Tau
figure();
hold on;

plot(uVals, q_all{1}, 'LineWidth', 2);
plot(uVals, q_all{2}, 'LineWidth', 2);

% plot(dataHard.displacement.value,dataHard.q.maxValue,'o','MarkerSize',10,'MarkerIndices',1:30:length(dataHard.q.maxValue),Color='#0072BD')
% plot(dataSoft.displacement.value,dataSoft.q.maxValue,'square','MarkerSize',10,'MarkerIndices',1:20:length(dataSoft.q.maxValue),Color='#D95319')

xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Hardening law (q)', 'Interpreter', 'latex');
% title('Exponential hardening law as a function of the damage evolution parameter');
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')




% Plot Force vs Displacement
figure();
hold on;
plot(uVals, F_all{1}, 'LineWidth', 2);
plot(uVals, F_all{2}, 'LineWidth', 2);

% plot(dataHard.displacement.value,dataHard.reaction,'o','MarkerSize',10,'MarkerIndices',1:30:length(dataHard.q.maxValue),Color='#0072BD')
% plot(dataSoft.displacement.value,dataSoft.reaction,'square','MarkerSize',10,'MarkerIndices',1:20:length(dataSoft.q.maxValue),Color='#D95319')

xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
% title('Reaction force as a function of the damage evolution parameter');
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')

% % Plot Displacement vs r
% fontsize(gcf,15,'points')
% figure();
% hold on;
% for i = 1:2
%     plot(uVals,max(tau_all{i},r0), 'LineWidth', 2);
% end
% xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
% ylabel('Damage evolution parameter (r)', 'Interpreter', 'latex');
% % title('Damage evolution parameter as a function of displacement');
% legend(labels, 'Location', 'best');
% grid on;

%%
clear; close all;

load('Swan/Resultats_definitius/1Elem-PF/1Elem_AT1.mat');
dataAT1_1elem = outputData;
load('Swan/Resultats_definitius/SENtraction-PF/AT1traction.mat');
dataAT1_Traction = outputData;
load('Swan/Resultats_definitius/SENshear-PF/AT1shear.mat');
dataAT1_Shear = outputData;
%------------------------------------------
load('Swan/Resultats_definitius/1Elem-PF/1Elem_AT2.mat');
dataAT2_1elem = outputData;
load('Swan/Resultats_definitius/SENtraction-PF/AT2traction.mat');
dataAT2_Traction = outputData;
load('Swan/Resultats_definitius/SENshear-PF/AT2shear.mat');
dataAT2_Shear = outputData;
%-------------------------------------------
load('Swan/Resultats_definitius/1Elem-CD/1ELM_LINEAR_SOFTENING_H_05_R0_10_R1_25.mat');
dataLinSoft_1elem = data;
load('Swan/Resultats_definitius/SEN Traction CD/LinearTraction.mat');
dataLinSoft_Traction = data;
load('Swan/Resultats_definitius/SEN Shear CD/LinearShear.mat');
dataLinSoft_Shear = data;
%-------------------------------------------
load('Swan/Resultats_definitius/1Elem-CD/1ELEM_TRACTION_EXP_SOFT.mat');
dataExpSoft_1elem = data;
load('Swan/Resultats_definitius/SEN Traction CD/SEN_TRACTION_EXP_SOFT_A_01_Qinf_2_R0_10.mat');
dataExpSoft_Traction = data;
load('Swan/Resultats_definitius/SEN Shear CD/ExpShear.mat');
dataExpSoft_Shear = data;
%-------------------------------------------
load('Swan/Resultats_definitius/1Elem-CD/1ELM_LINEAR_HARDENING_H_05_R0_10_R1_25.mat');
dataLinHard_1elem = data;
load('Swan/Resultats_definitius/SEN Traction CD/SEN_TRACTION_LINEAR_HARD_H_05_R0_10_R1_20.mat');
dataLinHard_Traction = data;
load('Swan/Resultats_definitius/SEN Shear CD/SEN_SHEAR_LINEAR_HARD_H_05_R0_10_R1_20.mat');
dataLinHard_Shear = data;
%-------------------------------------------
load('Swan/Resultats_definitius/1Elem-PF/1Elem-CD-AT1.mat');
dataCdAT1_1elem = data;
% dataCdAT1_Traction = data;
% dataCdAT1_Shear = data;
%-------------------------------------------
load('Swan/Resultats_definitius/1Elem-PF/1Elem-CD-AT2.mat');
dataCdAT2_1elem = data;
load('Swan/Resultats_definitius/SEN Shear CD/SEN_TRACTION_Q_AT2_W1_3750.mat');
dataCdAT2_Traction = data;
load('Swan/Resultats_definitius/SENshear-PF/SEN_SHEAR_Q_AT2_W1_3750.mat')
dataCdAT2_Shear = data;
%-------------------------------------------
% dataPfLinear_1elem = outputData;
% dataPfLinear_Traction = outputData;
% dataPfLinear_Shear = outputData;

%% 1 ELEM CD

figure ();
hold on

labels = {'Exponential Softening','Linear Softening', 'Linear Hardening'};
data_cd = [dataExpSoft_1elem,dataLinSoft_1elem,dataLinHard_1elem];

for i = 1:length(data_cd)
    plot(data_cd(i).displacement.value,data_cd(i).reaction,'LineWidth', 2);
end
xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
title('1Elem CD')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')

%----------------------------------------------------------------------------------------------
figure ();
hold on

labels = {'Exponential Softening','Linear Softening', 'Linear Hardening'};
data_cd = [dataExpSoft_1elem,dataLinSoft_1elem,dataLinHard_1elem];

for i = 1:length(data_cd)
    plot(data_cd(i).displacement.value,data_cd(i).damage.maxValue,'LineWidth', 2);
end
xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Damage ($d$)', 'Interpreter', 'latex');
title('1Elem CD')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')
%----------------------------------------------------------------------------------------------
%% 1 ELM PF

figure ();
hold on

labels = {'AT1','AT2','CD-AT1','CD-AT2'};
colorLabel = {'#0072BD' ,'#D95319'};
markerLabel = {'o','square'};

data_1elem_Pf_cd = [dataCdAT1_1elem,dataCdAT2_1elem];
data_1elem_Pf = [dataAT1_1elem,dataAT2_1elem];
for i = 1:length(data_1elem_Pf)
    plot(data_1elem_Pf(i).displacement.value,data_1elem_Pf(i).force,'LineWidth', 2);
end
for i = 1:length(data_1elem_Pf_cd)
        plot(data_1elem_Pf_cd(i).displacement.value, data_1elem_Pf_cd(i).reaction, ...
         'LineStyle', 'none', ...
         'Color', colorLabel{i}, ...
         'Marker', markerLabel{i}, ...
         'MarkerSize', 5, ...
         'MarkerIndices', 1:10:length(data_1elem_Pf_cd(i).reaction));
end

xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
title('1 Element PF')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')
%-------------------------------------------------
figure ();
hold on

labels = {'AT1','AT2','CD-AT1','CD-AT2'};
colorLabel = {'#0072BD' ,'#D95319'};
markerLabel = {'o','square'};

data_1elem_Pf = [dataAT1_1elem,dataAT2_1elem];
data_1elem_Pf_cd = [dataCdAT1_1elem,dataCdAT2_1elem];

for i = 1:length(data_1elem_Pf)
    plot(data_1elem_Pf(i).displacement.value,data_1elem_Pf(i).damage.maxValue,'LineWidth', 2);
end
for i = 1:length(data_1elem_Pf_cd)
        plot(data_1elem_Pf_cd(i).displacement.value, data_1elem_Pf_cd(i).damage.maxValue, ...
         'LineStyle', 'none', ...
         'Color', colorLabel{i}, ...
         'Marker', markerLabel{i}, ...
         'MarkerSize', 5, ...
         'MarkerIndices', 1:10:length(data_1elem_Pf_cd(i).reaction));
end

xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Damage ($d$)', 'Interpreter', 'latex');
title('1 Element PF')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')

%% Traction CD

figure ();
hold on

labels = {'Exponential Softening','Linear Softening', 'Linear Hardening'};
data_cd = [dataExpSoft_Traction,dataLinSoft_Traction,dataLinHard_Traction];

for i = 1:length(data_cd)
    plot(data_cd(i).displacement.value,data_cd(i).reaction,'LineWidth', 2);
end
xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
title('SEN Traction CD')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')
hold off
%----------------------------------------------------------------------------------------------
% for i = 1:length(data_cd)
%     figure();
%     hold on
%     data_cd(i).damage.field.plot;
%     title (labels{i})
%     hold off;
% end

%% Traction PF

figure ();
hold on

labels = {'AT1','AT2','CD-AT2'};
colorLabel = {'#D95319'};
markerLabel = {'square'};

data_1elem_Pf = [dataAT1_Traction,dataAT2_Traction];
data_1elem_Pf_cd = [dataCdAT2_Traction];

for i = 1:length(data_1elem_Pf)
    plot(data_1elem_Pf(i).displacement.value,data_1elem_Pf(i).force,'LineWidth', 2);
end
for i = 1:length(data_1elem_Pf_cd)
        plot(data_1elem_Pf_cd(i).displacement.value, data_1elem_Pf_cd(i).reaction, ...
         'LineStyle', 'none', ...
         'Color', colorLabel{i}, ...
         'Marker', markerLabel{i}, ...
         'MarkerSize', 5, ...
         'MarkerIndices', 1:10:length(data_1elem_Pf_cd(i).reaction));
end

xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
title('SEN Traction PF')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')

%% SHEAR CDM

figure ();
hold on

labels = {'Exponential Softening','Linear Softening', 'Linear Hardening'};
data_cd = [dataExpSoft_Shear,dataLinSoft_Shear,dataLinHard_Shear];

for i = 1:length(data_cd)
    plot(data_cd(i).displacement.value,data_cd(i).reaction,'LineWidth', 2);
end
xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
title('SEN Shear CD')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')
hold off
%----------------------------------------------------------------------------------------------
% for i = 1:length(data_cd)
%     figure();
%     hold on
%     data_cd(i).damage.field.plot;
%     title (labels{i})
%     hold off;
% end

%% Shear PF

figure ();
hold on

labels = {'AT1','AT2','CD-AT2'};
colorLabel = {'#D95319'};
markerLabel = {'square'};

data_1elem_Pf = [dataAT1_Shear,dataAT2_Shear];
data_1elem_Pf_cd = [dataCdAT2_Shear];

for i = 1:length(data_1elem_Pf)
    plot(data_1elem_Pf(i).displacement.value,data_1elem_Pf(i).force,'LineWidth', 2);
end
for i = 1:length(data_1elem_Pf_cd)
        plot(data_1elem_Pf_cd(i).displacement.value, data_1elem_Pf_cd(i).reaction, ...
         'LineStyle', 'none', ...
         'Color', colorLabel{i}, ...
         'Marker', markerLabel{i}, ...
         'MarkerSize', 5, ...
         'MarkerIndices', 1:10:length(data_1elem_Pf_cd(i).reaction));
end

xlabel('Displacement ($\varepsilon$)', 'Interpreter', 'latex');
ylabel('Reaction ($\sigma$)', 'Interpreter', 'latex');
title('SEN Shear PF')
legend(labels, 'Location', 'best');
grid on;
fontsize(gcf,15,'points')
