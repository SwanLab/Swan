% Computes one V-cycle
function [v] = MG(nx_current, ny_current, A_h, f_h, v_h)

if nx_current == 5 % This is the number of cells in x in coarsest mesh
    v = A_h\f_h;
else
nx_coarse = (nx_current+1)/2; % Number of x nodes in next coarsest mesh
ny_coarse = (ny_current+1)/2; % Number of y nodes in next coarsest mesh
I_2h_h = construct_I(nx_coarse, ny_coarse); % Interpolation matrix
R_h_2h = construct_R(I_2h_h); % Restriction matrix
A_2h = R_h_2h * A_h * I_2h_h; % A_2h = R*A_h*I

v_h = gauss_seidel_fixed(A_h, f_h, v_h, 3); % 3 Gauss seidel steps

r_2h = R_h_2h * (f_h - A_h*v_h); % Projection of residual to coarse mesh

% Recursive implementation of calculating error
e_2h = MG(nx_coarse, ny_coarse, A_2h, r_2h, zeros(size(r_2h))); 

% Projecting error from coarse to fine mesh
e_h = I_2h_h*e_2h;

% Initial guess correction to coarse mesh
v = v_h + e_h;
end

    