% MATLAB Script to compute the Loading Function L
clear; clc;

%% 1. User Inputs
tau3_o   = 15;     % Interlaminar strength in Mode I (so 3)
tau_sh_o = 30;     % Interlaminar strength in shear mode (sosh)
G_Ic     = 260;    % Pure mode I fracture toughness (GIc)
G_IIc    = 1002;    % Pure mode II/shear fracture toughness (GIIc)
eta      = 2;    % Mixed-mode interaction parameter (g)

% New Explicit State Inputs
lambda   = 0.05;   % Current mixed-mode displacement jump norm
Delta_sh = 0.03;   % Displacement jump in shear mode

% Additional required simulation parameter
k        = 1e5;    % Interface penalty stiffness

%% 2. Calculate Mixed-Mode Ratio (Equation 11)
if lambda == 0
    B = 0; % Prevent division by zero at zero loading
else
    B = (Delta_sh^2) / (lambda^2);
end

% Ensure B stays physically bounded between 0 and 1 due to numerical precision
B = max(0, min(1, B)); 

%% 3. Intermediate Calculations (Equations 14 & 15)
% Equation (14): Mixed-mode interface strength (tau_o)
tau_o = sqrt( (tau3_o)^2 + ( (tau_sh_o)^2 - (tau3_o)^2 ) * B^eta );

% Equation (15): Mixed-mode fracture toughness (G_c)
G_c = G_Ic + (G_IIc - G_Ic) * B^eta;

%% 4. Calculate the Loading Function (Equation 16)
numerator   = 2 * G_c * (k * lambda - tau_o);
denominator = lambda * (2 * k * G_c - (tau_o)^2);

% Handle edge case where lambda or denominator is zero
if denominator == 0 || lambda == 0
    L_term = 0; 
else
    L_term = numerator / denominator;
end

% Equation (16): L = min { Term , 1 }
L = min(L_term, 1);

%% 5. Display the Output
fprintf('--- Results ---\n');
fprintf('Calculated Mixed-mode ratio (B): %.4f\n', B);
fprintf('Mixed-mode strength (tau_o):     %.4f\n', tau_o);
fprintf('Mixed-mode toughness (G_c):      %.4f\n', G_c);
fprintf('Loading Function (L):            %.4f\n', L);