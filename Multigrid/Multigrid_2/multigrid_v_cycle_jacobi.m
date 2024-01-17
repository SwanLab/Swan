
function u = multigrid_v_cycle_jacobi(u, f, L, Ncycles, Npre, Npost)
% Performs Ncycles V-cycles of multigrid with bilinear interpolation and Jacobi smoother

% Define relaxation function
omega = 2/3; % Jacobi relaxation parameter
%relax = @(u, f, h) jacobi(u, f, h, omega, Npre, Npost);
relax = @(u,f,h) jacobi_2d(u, f, h, Npre, Npost);

% Define restriction function
restrict = @(u) restriction(u);


% Define interpolation function
interp = @(u) bilinear_interpolation(u);


% Define residual function
res = @(u, f, h) residual(u, f, h);

% Define exact solution function (for testing)
exact = @(x, y) exact_u(x, y);

% Define initial grid spacing and number of levels
h = 1/(L-1);
num_levels = log2(L-1);
%disp(num_levels)

% define rhs
%f = @(x,y)   rhs_func(x,y);

% Define hierarchy of grids
U = cell(num_levels, 1); % stores the numerical solutions
F = cell(num_levels, 1); % stores the right-hand sides
R = cell(num_levels, 1); % stores the residuals
U{1} = u;
F{1} = f;
for k = 2:num_levels
    n = (L-1)/2^(k-1) + 1; % number of grid points in each direction
    %disp(n)
    U{k} = zeros(n, n); % initialize numerical solution
    F{k} = zeros(n, n); % initialize right-hand side
end

% Perform V-cycles
for cycle = 1:Ncycles
    % Downward sweep
    for k = 1:num_levels-1
        % Pre-smoothing
        U{k} = relax(U{k}, F{k}, h);
       
        % Compute residual and restrict to coarser grid
        R{k} = res(U{k}, F{k}, h);
        %disp(R{k})
        %disp(size(R{k}))
        F{k+1} = restrict(R{k});
        %disp(F{k+1})
       
        % Initialize solution guess for next level
        U{k+1} = zeros(size(F{k+1}));
    end
   
    % Bottom-level correction
        U{num_levels} = direct_solve(F{num_levels}, h);
        %disp(U{num_levels})
  
    % Upward sweep
    for k = num_levels-1:-1:1
        % Interpolate and add correction to solution
        E = interp(U{k+1});
        %disp(E)

        U{k} = U{k} + E;
       % disp(U{k});
        
        % Post-smoothing
        U{k} = relax(U{k}, F{k}, h);
        
        
        %M = cell2mat(U(1));
        %disp(M)
    
    end
    % Compute residual and error
    R{1} = res(U{1}, F{1}, h);
    err = norm(R{1}(:), 2)/ sqrt((L-1)^2);
    fprintf('Cycle %d: error = %e\n', cycle, err);
    
   
end
    L = 129; % number of grid points in each direction
    [X, Y] = meshgrid(linspace(0, 1, L), linspace(0, 1, L));
    
    %f = rhs_func(X,Y);


    figure;
    %subplot(1, 2, 1);
    surf(X, Y, U{k});
    title('Numerical solution');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    %subplot(1, 2, 2);
 
end
    

