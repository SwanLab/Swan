  clc;
  clear;
  close;
%% 2D heat equation with forcing  
% Solve for $u(x,y,t)$ such that
%
% $$u_{t}-\Delta u= (\pi\sin(2\pi t)-\pi^{2}\cos(2\pi t))\sin(\pi x)\sin(\pi y),\textnormal{on }\Omega,t>0,$$
%
% $$u=0,\textnormal{on }\partial\Omega,t>0,$$
%
% $$u(x,y,t)=-\frac{1}{2}\sin(\pi x)\sin(\pi y)$$
%
% where $\Omega={(x,y):0<x<2,0<y<1}$ and $\partial\Omega$ is boundary of $\Omega$
% 
% # Rectangular domain: [x_start, x_end] X [y_start, y_end]
% # Rectangular mesh with spacing h in x and y directions (ndim = 2)
% # Mesh connectivity stored in inode
% # Tranformation from x,y directions to zi,eta directions
% # nnodes = 4 nodes in rectangular element
% # Dirichlet condition with 0 at boundary (id = 0 at boundary, 1 elsewhere)
% # Gaussian quadrature integration with num_gauss_quad_points quadrature points
% # Quadrature points stored as (zi_g, eta_g), 2D array of size (num_gauss_quad_points, ndim)
% # Shape functions in zi, eta system defined in compute_shape.m
% # Derivative of shape functions defined in compute_df

%% Set true solution for comparison later
% $$u(x,y,t)=-\frac{1}{2}\cos(2\pi t)\sin(\pi x)\sin(\pi y)$$
set_true_sol = @(t,x,y) -0.5*(cos(2*pi*t)*sin(pi*x).*sin(pi*y));

%% Set force
% $$f(x,y,t)=(\pi\sin(2\pi t)-\pi^{2}\cos(2\pi t))\sin(\pi x)\sin(\pi y)$$
set_force = @(t,x,y) (pi*sin(2*pi*t) - ...
    pi*pi*cos(2*pi*t))*(sin(pi*x)*sin(pi*y));

%% Set initial condition
% $$u(x,y,t)=-\frac{1}{2}\sin(\pi x)\sin(\pi y)$$           
set_IC = @(x,y) -0.5*sin(pi*x).*sin(pi*y);

%% Input parameters
xstart = 0.0; % x0
xend = 2.0; % x1
ystart = 0.0; % y0
yend = 1.0; % y1
h = 0.25; % Spacing
tend = 4.0; % Final time
theta=0.5; % 0 <= theta <= 1 | 1: Explicit; 0: Implcit
dt = 0.1; % Time step
sol_method = 1; % 1 = FEM; 2 = FDM (5 point stencil)
solver = 3; % 0: MATLAB default solver
            % 1: Conjugate gradient
            % 2: Gauss_seidel 
            % 3: Multigrid with Gauss seidel iteration upto ...
            % number of nodes in x direction become 5

%% Finite element parameters
nnode = 4; % Number of nodes in element
ndim = 2; % 2 directions x & y

%% Parameters for numerical integration - Setting quadrature points
num_quad_quad_points = 4; % Number of gauss points
if (sol_method == 1) % If method is FEM
    zi_1 = [-1.0/sqrt(3), 1.0/sqrt(3)]; % z direction gauss points
    eta_1 = [-1.0/sqrt(3), 1.0/sqrt(3)]; % eta direction gauss points
    w_1 = ones(num_quad_quad_points,1); % Weights for integration for zi direction
    w_2 = ones(num_quad_quad_points,1); % Weights for integration for eta direction
elseif(sol_method == 2) % If method is 5 point FDM
    zi_1 = [-1.0, 1.0]; % z direction gauss points
    eta_1 = [-1.0, 1.0]; % eta direction gauss points
    w_1 = 0.25*ones(num_quad_quad_points,1); % Weights for integration for zi direction
    w_2 = 0.25*ones(num_quad_quad_points,1); % Weights for integration for eta direction
end
[zi_g, eta_g] = meshgrid(zi_1,eta_1);
zi_g = zi_g(:);
eta_g = eta_g(:);

%% Set the mesh: These variables will be used later
num_nodes_x = fix((xend-xstart)/h)+1; % Nodes in x direction
num_nodes_y = fix((yend-ystart)/h)+1; % Nodes in y direction
nj = num_nodes_x*num_nodes_y; % Total number of nodes
nelem = (num_nodes_x-1)*(num_nodes_y-1); % Total number of nodes

x_points = linspace(xstart,xend,num_nodes_x);
y_points = linspace(ystart,yend,num_nodes_y);
[y_mesh,x_mesh] = meshgrid(y_points,x_points);
x = x_mesh(:);
y = y_mesh(:);

y_old = set_IC(x,y); %IC

%% Set boundary condition using array id
id = zeros(num_nodes_x,num_nodes_y);
id(1,:) = 1;
id(end,:) = 1;
id(:,1) = 1;
id(:,end) = 1;
id = id(:);

%% Set mesh connectivity in array inode
inode = zeros(nnode,nelem);
iter1 = 1;
iter2 = 0;
for element_num=1:nelem
    inode(1,element_num) = iter1 + iter2*(num_nodes_x);
    inode(2,element_num) = inode(1,element_num) + 1;
    inode(3,element_num) = inode(2,element_num) + num_nodes_x;
    inode(4,element_num) = inode(3,element_num) - 1;
    iter1 = iter1 + 1;
    if(mod(element_num,num_nodes_x-1) == 0)
        iter1 = 1;
        iter2 = iter2 + 1;
    end
end

%% Time integration: Loop over every element
t0 = 0.0;
t1=0.0;
iter_array = zeros(fix(tend/dt),1);
iter_num = 1;
while(abs(t1 - tend) > dt/4.0)
    t1 = t0+dt;
    global_k = zeros(nj,nj);
    global_m = zeros(nj,nj);
    global_f = zeros(nj,1);
    
    for element_num=1:nelem
        
        xt = zeros(ndim,nnode); % x and y points of an element
        % Obtain the x and y coordinates from global connectivity
        xt(1,:) = x(inode(:,element_num));
        xt(2,:) = y(inode(:,element_num));
        
        ke = zeros(nnode,nnode); % local element stiffness matrix
        me = zeros(nnode,nnode); % local element mass matrix
        fe = zeros(nnode,1); % Local force vector
        % Integrate using gauss points
        for gauss_points=1:num_quad_quad_points
            grad_N = compute_dfx(xt, zi_g(gauss_points), ...
                eta_g(gauss_points), nnode, ndim);
            N = compute_shape(zi_g(gauss_points), eta_g(gauss_points), ...
                nnode);
            
            J_dx_dy = det(transpose(compute_dj(xt,zi_g(gauss_points), ...
                eta_g(gauss_points), nnode, ndim)));
            
            Nx = grad_N(:,1);
            Ny = grad_N(:,2);
            
            ke = ke + (Nx*Nx' + Ny*Ny')*  w_1(gauss_points) * ...
                w_2(gauss_points) * J_dx_dy; % local stiffness ke
            me = me + (N*N')*  w_1(gauss_points) * ...
                w_2(gauss_points) *J_dx_dy;
            
              fe = fe + (theta*N*set_force(t0, xt(1,:)*N, xt(2,:)*N) + ...
                  (1-theta)*N*set_force(t1, xt(1,:)*N, xt(2,:)*N)) *  ...
                  w_1(gauss_points) * w_2(gauss_points) *J_dx_dy;
        end
        
        % Assembly
        IG_rows = inode(:,element_num);
        for il = 1:nnode
            IG = IG_rows(il);
            for jl = 1:nnode
                JG = IG_rows(jl);
                global_k(IG,JG) = global_k(IG,JG) + ke(il,jl); % Global stiffness K
                global_m(IG,JG) = global_m(IG,JG) + me(il,jl); % Global mass M
            end
            global_f(IG) = global_f(IG) + fe(il); % Global force F
        end
    end
    
    A = global_m + dt*(1-theta)*global_k;
    b = global_f*dt + (global_m - dt*theta*global_k)*y_old;
    
    % Apply BC
    % Dirichlet condition = 0 implies no terms go to RHS
    for eq_num=1:nj
        if (id(eq_num) ~= 0)
            A(:,eq_num) = 0.0;
            A(eq_num,:) = 0.0;
            A(eq_num,eq_num) = 1.0;
            b(eq_num) = 0.0;
        end
    end
    
    A = sparse(A);
    
    if (solver == 0) % Default MATLAB solver
        y_new = A\b;
    elseif(solver == 1) % Conjugate gradient
        [y_new, cg_iter] = cg(A,b);
        iter_array(iter_num,1) = cg_iter;
        iter_num = iter_num + 1;
    elseif(solver == 2) % Gauss seidel
        [y_new, gs_iter] = gauss_seidel(A,b);
        iter_array(iter_num,1) = gs_iter;
        iter_num = iter_num + 1;
    elseif(solver == 3) % Multigrid
        iter_mg = 0;
        [y_new, mg_iter] = MG_full_cycle(num_nodes_x, num_nodes_y, A, b);
        iter_array(iter_num,1) = mg_iter;
        iter_num = iter_num + 1;
    end
    
    y_old = y_new;
    
    t0 = t0 + dt
end
solut_mesh = reshape(y_old,[num_nodes_x, num_nodes_y]);

%% Post-processing
subplot(2,2,1)
surf(x_mesh,y_mesh,solut_mesh)
xlabel('X');
ylabel('Y');
str = sprintf('Numerical solution: h = %4.3f (%d elements), dt = %4.3f ', ...
    h, nelem, dt);
title(str)
axis equal;
colorbar;
caxis([-0.5, 0.5]);
colorbar('off')

true_solut_mesh = set_true_sol(tend, x_mesh, y_mesh);
subplot(2,2,2)
surf(x_mesh,y_mesh,true_solut_mesh)
xlabel('X');
ylabel('Y');
title('True solution')
axis equal;
colorbar('Ticks', linspace(-0.5,.5,11));
caxis([-0.5, 0.5]);

error_mesh = true_solut_mesh - solut_mesh;
subplot(2,2,3);
surf(x_mesh,y_mesh,error_mesh)
xlabel('X');
ylabel('Y');
title('Error')
axis equal;
c = colorbar;

subplot(2,2,4);
plot(iter_array, 'ro');
ytickformat('%,.0f');
if solver == 1
    str = sprintf('# CG iterations: h = %4.3f (%d elements), dt = %4.3f ', ...
    h, nelem, dt);
elseif solver == 2
    str = sprintf('# GS iterations: h = %4.3f (%d elements), dt = %4.3f ', ...
    h, nelem, dt);
elseif solver == 3
    str = sprintf('# V-cycles (1 cycle has 3 GS iterations): h = %4.3f (%d elements), dt = %4.3f ', ...
    h, nelem, dt);
end
title(str);
yticks(1:max(iter_array))
xlabel('Timestep')
ylabel('Iterations')

if (sol_method == 1)
    subtitle((['FEM scheme \theta=',num2str(theta), ...
        ' (\theta=0: Implicit \theta=1: Explicit): Comparison with true solution']));
elseif (sol_method == 2)
    suptitle((['(Multigrid) FDM (5 point stencil) \theta = ',num2str(theta),': Comparison with true solution']));
end
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 26 20])
set(gcf,'PaperOrientation','landscape');
print(gcf, '-dpdf', 'dummy.pdf')
