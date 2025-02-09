%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT :                                                            |
%| PROPERTY  :                                                            |
%| LICENCE   :                                                            |
%| CONTACT   :                                                            |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtStkTranslatingSphere.m                     |
%|    #    |   VERSION    : 0.10                                          |
%|   _#_   |   AUTHOR(S)  : Luca Berti                                    |
%|  ( # )  |   CREATION   : 25.12.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+
% TEST CASE - TRANSLATING SPHERE IN THE Z DIRECTION
% 
% Testcase to verify the BEM Stokes solver: in presence of a sphere, we
% want to make the solution of it to radiate. We simulate a translating 
% sphere with given velocity U and radius a. 
% 
% We proceed as follows: 
% 1) we compute the solution of the single-layer equation U = G * stress,
% where G = Stokeslet + n \otimes n and "stress" is the unknown
% 2) once we know the stresses on the sphere surface, we evaluate the
% solution coming from Happel-Brenner, "Low Reynolds number hydrodynamics"
% (pag. 119 and formulas 4-17.17 and 4-17.18) transformed into Cartesian 
% coordinates (p.507 A-15.29 for the transformation) for comparison 
% with the numerical solution we get from U_num = -G_eval*stresses

clear
close all
clc

% Gypsilab path
run('../../../addpathGypsilab.m')

% PHYSICAL PARAMETERS
W = 2; % Translational velocity of the sphere - z component
viscosity = 1; % Fluid viscosity
sphere_translational_velocity =[0;0;W]; % U
sphere_radius = 1;
alpha = 0; % coefficient weighting the normal \otimes normal correction on the single-layer operator matrix
pmin = 6;
pmax = 9; % total number of iterations for the convergence test
mesh_size = zeros(pmax,1);
L2error = zeros(pmax,1);
Linferror = zeros(pmax,1);


tol = []; % CHOOSE WHETHER TO USE H-MATRICES OR NOT
            % [] -> no H-matrix
            % 1e-5 (example) -> H-matrix with prescribed tolerance
test_local_refinement = false; % CHOOSE WHETHER TO TEST A NON-LOCALLY
                               % REFINED MESH (true or false)
if isempty(tol) == 1
    disp('Using BEM without H-matrix');
else 
    disp('Using H-matrices')
end
if test_local_refinement == false
    disp('Not testing local refinement')
else
    disp('Testing local refinement')
end
    

%% CONVERGENCE TEST - possibility to use H-matrices and locally refined meshes

for p = pmin:pmax
    % CONSTRUCTION OF THE COMPUTATIONAL DOMAIN AND QUANTITIES
    %tic
    2^p
    mesh = mshSphere(2^p,sphere_radius);
    if test_local_refinement==true
        mesh = refine(mesh,@(z) 50);
        % Project the newly computed points to the sphere surface
        for i = 1:length(mesh.vtx)
            mesh.vtx(i,:) = sphere_radius*mesh.vtx(i,:)./(norm(mesh.vtx(i,:)));
        end
    end
    Gamma = dom(mesh,3);
    phi = fem(mesh,'P1');
    N_unknown = size(phi.unk,1);
    
    % CONSTRUCTION OF THE EVALUATION POINTS
    N_eval = 1000;
    mesh_eval = mshSphere(N_eval,2*sphere_radius);

    % Compute the single-layer G and mass M matrices
    tic
    [G,M] = stokesletEQ(Gamma,phi,viscosity,alpha,tol);
    toc
    % Compute the right-hand side of the single-layer problem
    U = rhs_computation(M,sphere_translational_velocity);
    
    % Compute the boundary stresses (either with GMRES either with a direct
    % method)
    %stress_sol = gmres(-G,U,[],1e-6,100);
    stress_sol = -G\U;
    
    % Compute the Stokeslet radiation matrix and the numerical solution U_sol
    % in correspondence of the evaluation points
    G_eval = StokesletRAD(Gamma,phi,viscosity,mesh_eval);
    U_sol = (-G_eval*stress_sol);
    
    % Compute the exact solution of the translating sphere problem from the
    % formulas given in the reference
    U_exact = exact_translating_sphere(sphere_translational_velocity,sphere_radius,mesh_eval);
    
    % L2 and L\infty norm of the error
    disp(["L2 norm of the error is ", norm(U_exact-U_sol)/norm(U_exact)])
    disp(["L \infty norm of the error is ", norm(U_exact-U_sol,'inf')/norm(U_exact,'inf')])
    
    temp = mesh.stp;
    mesh_size(p) = temp(2);
    L2error(p) = norm(U_exact-U_sol)/norm(U_exact);
    Linferror(p) = norm(U_exact-U_sol,'inf')/norm(U_exact,'inf');
    %toc
end

% % Computation of the Stokes resistance force
% stokes_law = -6*pi*viscosity*sphere_radius*W;
% numerical_law = sum(stress_sol(2*N_unknown+1:end,1)*4*pi*sphere_radius^2/N_unknown);
% error_stokeslaw = abs(numerical_law-stokes_law)/abs(stokes_law)

% Experimental order of convergence
for p=pmin+1:pmax
    eoc(p) = log(L2error(p)/L2error(p-1))/log(mesh_size(p)/mesh_size(p-1));
end

% Plot of convergence graph
loglog(mesh_size,L2error,'b+-',mesh_size,mesh_size,'k',mesh_size,(mesh_size*10^-0.5).^2,'k--')
legend('L2error','slope 1','slope 2')
title('Log-log plot of mesh size vs. L^2 error')
xlabel('mesh size')
ylabel('L^2 error')





disp('~~> Michto gypsilab !')



