% Copyright (c) 2021 STIJN KOPPEN
%
% function flexure( nelx, nely, doc, dof, emax)
%
% Author1:  Stijn Koppen (s.koppen@tudelft.nl)
% Author2:  Matthijs Langelaar
% Author3:  Fred van Keulen
% Date:     1 June 2021
% Updated:  2 November 2022
%
% Function:     flexure
% Description:  Synthesis of a 2D flexure using topology optimization
% Parameters:   nelx    int                     number of elements in x-direction
%               nely    int                     number of elements in y-direction
%               doc     string/list of strings  degrees of constraint (DOCs)
%               ("tx", "ty" and/or "rz")
%               dof     string/list of strings  degrees of freedom (DOFs)
%               ("tx", "ty" and/or "rz")
%               emax    float/list of floats    maximum strain energy for
%               DOFs
%
% Description of degrees:   "tx": translation in x-direction (horizontal)
%                           "ty": translation in y-direction (vertical)
%                           "rz": rotation about z-direction (point of
%                           rotation in middle)
%
% Examples of usage:
% To synthesize a 2D flexure with 100 x 100 finite elements, stiff in y-translation, while compliant (strain energy <= 1 Nm) in
% x-translation:    flexure( 100, 100, "ty", "tx", 1.0)
%                   flexure( 100, 100, ["ty"], ["tx"], 1.0)
%
% To synthesize a 2D flexure, stiff in x- and y-translation, while compliant (strain energy <= 0.1 Nm) 
% about the z-axis:    flexure( 100, 100, ["tx", "ty"], "rz", 0.1)
%
% To synthesize a 2D flexure, stiff in x-translation, while compliant (strain energy <= 1 Nm) 
% in y-translation AND about the z-axis (strain energy <= 0.1 Nm):    flexure( 100, 100, "tx", ["ty", "rz"], [1, 0.1])
%
% Note: maximum strain energy can easily be calculated from stiffness and vice
% versa.
% For translational degrees E = 0.5 * F * u = 0.5 * k * u^2 (Nm), with force F and linear stiffness k = F/u (N/m). In this code a unit
% displacement is applied (u = 1 m), so E = 0.5 * k.
% For rotational degree E = 0.5 * T * theta = 0.5 * k * theta^2 (Nm), with torque T (Nm) and rotational
% stiffness k = T/theta (Nm). In this code a unit rotation is applied
% (theta = 1), so E = 0.5 * k.
%
% Finite element node and element numbering for nelx = 3 and nely = 2 (6 elements, 12 nodes):
%
%    o - x
%    |
%    y  1 - 4 - 7 - 10
%       | 1 | 3 | 5 |
%       2 - 5 - 8 - 11
%       | 2 | 4 | 6 |
%       3 - 6 - 9 - 12
%%%% A TOPOLOGY OPTIMIZATION CODE FOR FLEXURE SYNTHESIS (2021) %%%%
function flexure(nelx, nely, doc, dof, emax)
    %% PREPROCESSING
    ldoc = length(doc); % Number of DOCs
    ldof = length(dof); % Number of DOFs
    
    % mdoc and mdof are used to compare input (doc, dof) to degrees
    mdoc = zeros(ldoc,3);
    mdof = zeros(ldof,3);
    
    degrees = ["tx","ty","rz"];
    for i=1:ldoc; mdoc(i,:) = strcmp(degrees,doc(i)); end % Compare DOCs to degrees
    for i=1:ldof; mdof(i,:) = strcmp(degrees,dof(i)); end % Compare DOFs to degrees
    mdoc = max(mdoc,[],1) == 1;
    mdof = max(mdof,[],1) == 1;
    
    deg = find(max([mdoc; mdof])); % Degrees that are in doc or dof
    ndeg = length(deg); % Number of active degrees
    
    assert(all((mdoc+mdof) < 2),'overlap between DOC and DOF');
    assert(ldoc >= 1 & ldoc <= 2,'set of DOC too small/big');
    assert(ldof >= 1 & ldof <= 2,'set DOF too small/big');
    assert(ndeg >= 2 & ndeg <=3,'incorect no of rhs');
    
    %% CONSTANTS
    E = 1000; % Youngs' Modulus (Pa)
    dz = 1; % Element thickness (m)
    nu = 0.3; % Poisson ratio
    eps = 1e-9; % Relative Young's Modulus: Evoid/Esolid
    rmin = 2.0; % Density filter radius ( expressed in no. of elements )
    penal = 3.0; % Penalization factor for Modified SIMP interpolation function
    
    %% PREPARE MESH
    n = nelx*nely; % Number of variables, number of finite elements
    nody = nely + 1; % Number of nodes in y-direction
    nodx = nelx + 1; % Number of nodes in x-direction
    ndofy = 2*nody; % Number of nodal dof in y-direction
    ndof = ndofy*nodx; % Number of nodal dof
    ndpe = 8; % Number of nodal dof per element
    
    %% PREPARE FEA
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
    KE = E*dz/24*(1-nu^2)*([A11 A12; A12' A11] + nu*[B11 B12; B12' B11]);
    nodenrs = reshape(1:nodx*nody, nody, nodx);
    edofVec = reshape(2*nodenrs(1:end-1, 1:end-1) + 1, n, 1);
    edofMat = repmat(edofVec, 1, ndpe ) + repmat([0 1 2*nely + [2 3 0 1] -2 -1], n, 1);
    iK = reshape(kron(edofMat, ones(ndpe, 1))', n*ndpe^2, 1);
    jK = reshape(kron(edofMat, ones(1, ndpe))', n*ndpe^2, 1);
    
    %% PREPARE FILTER
    [dy,dx] = meshgrid(-ceil(rmin) + 1:ceil(rmin)-1, -ceil(rmin) + 1:ceil(rmin)-1);
    h = max(0, rmin - sqrt(dx.^2 + dy.^2));
    Hs = conv2(ones(nely, nelx), h, 'same');
    
    %% BOUNDARY CONDITIONS
    U = zeros(ndof, 3); % Initialization of displacement vectors
    
    % Axes: x, y, z with corresponding displacements u, v, w
    % ty: bottom nodes fixed, top nodes u=0, v=1
    % tx: bottom nodes fixed, top nodes u=1, v=0
    % rz: bottom nodes fixed, top nodes u=1, v=1-2*x/nelx, with x=0:nelx
    % over top nodes
    % NOTE: by changing the relative magnitude of prescribed displacements
    % for rz (u, v) one can change the point of rotation
    % NOTE: prescribed displacements are in meters, adjust according to
    % application
    
    U(2:ndofy:end, 2) = 1; % ty
    U(1:ndofy:end, [1 3]) = 1; % tx, rz
    U(2:ndofy:end, 3) = linspace(1, -1, nelx+1); % rz

    top = union(1:ndofy:ndof, 2:ndofy:ndof); % Get nodal dofs of top
    bottom = union(ndofy-1:ndofy:ndof, ndofy:ndofy:ndof); % Get nodal dofs of bottom
    fixed = union(top, bottom); % Fixed nodal dofs at bottom and top
    free = setdiff(1:ndof, fixed); % Free nodal dofs
    
    
    %% INITIALIZE OPTIMIZATION
    maxiter = 200;
    toldx = 5e-4; % Minimum mean variable change
    tolkkt = 5e-4; % KKT condition
    tolg = 1e-4; % Primal feasibility
    gscale = 0.01; % Constraint scaling
    % NOTE: ensure max constraint violation abs(g(x)) < 10
    %x = rand(n, 1); % Initial design (random)
    % Change to x = 0.2*ones(n,1) for homogeneous (feasible) initial design
    x = 0.5*ones(n,1);
    xold1 = x;
    xold2 = x;
    mlinit = 0.5; % Initial movelimit
    mlincr = 1.2; % Movelimit increase factor
    mldecr = 0.7; % Movelimit decrease factor
    movelimit = mlinit*ones(n, 1); % Initial move limit
    strain_energy = zeros(1, 3);
    dc = zeros(nely, nelx, 3); % Sensitivities
    
    %% SEQUENTIAL APPROXIMATE OPTIMIZATION LOOP
    iter = 0;
    while iter < maxiter
        iter = iter + 1;
        x = reshape(x, nely, nelx);
    
        %% FORCE SYMMETRY
        x = (x + flip(x,1))/2; % Symmetrize around x-axis
        x = (x + flip(x,2))/2; % Symmetrize around y-axis (elimination of round-off errors)
        
        %% DENSITY FILTER
        xPhys = conv2(reshape(x, nely, nelx), h, 'same')./Hs; % Filter design variables 
        % (using convolution function, see top71.m)
        
        %% MATERIAL INTERPOLATION
        E = eps + (1 - eps)*xPhys(:).^penal; % Material Interpolation (SIMP)
        
        %% STIFFNESS MATRIX ASSEMBLY
        sK = reshape(KE(:)*E', n*ndpe^2, 1); % Generate matrix entries
        K = sparse(iK, jK, sK); % Build sparse matrix
        K = (K + K')/2; % Ensure symmetry (elimination of round-off errors)
        
        %% SOLVE SYSTEM OF EQUATIONS
        U(free, deg) = K(free, free) \ (-K(free, fixed)*U(fixed, deg)); % x = A\b
        % NOTE: only solve for relevant degrees
        strain_energy(deg) = 0.5*sum(U(:, deg).*(K*U(:, deg))); % Calculate strain energies
        
        %% RESPONSES
        if iter == 1; alpha = 1./strain_energy; alpha(3) = alpha(3)*2;end % Scaling with respect to strain energy of first iteration
        f = -1/ldoc*sum( alpha(mdoc).*strain_energy(mdoc) ); % Objective
        g = gscale*(strain_energy(mdof)./emax - 1); % Constraint(s)
        
        % Alternative, one can scale the constraints wrt the first
        % iteration, thus emax = 1 -> g^{0} = 0
        % g = alpha(mdof).*strain_energy(mdof)./emax - 1;
        % Use 0.001 < emax < 100 (indication).
        % In that case, emax becomes a relative dimensionless number.
        
        %% SENSITIVITY ANALYSIS
        for i = 1:ndeg
            Ui = U(:,deg(i)); % Dummy variable Ui contains displacement field i
            ce = 0.5*reshape( sum( ( Ui(edofMat)*KE).*Ui(edofMat), 2 ), nely, nelx ); % Element strain energy
            dc(:,:,deg(i)) = penal*(1-eps)*xPhys.^(penal-1).*ce; % Sensitivity of strain energy i wrt filtered variables
            dc(:,:,deg(i)) = conv2(dc(:,:,deg(i))./Hs,h,'same'); % Sensitivity of strain energy i wrt symmetrized variables
            dc(:,:,deg(i)) = (dc(:,:,deg(i)) + flip(dc(:,:,deg(i)),2))/2;
            dc(:,:,deg(i)) = (dc(:,:,deg(i)) + flip(dc(:,:,deg(i)),1))/2;
            % NOTE the reversed order (first y-axis then x-axis)
        end        
        df = -1/ldoc*sum(alpha(mdoc).*reshape( dc(:,:,mdoc), n, ldoc),2); % Sensitivities of objective
        dg = gscale*(reshape( dc(:,:,mdof), n, ldof)./emax); % Sensitivities of constraints
        
        % Alternative use (see lines 170-173)
        % dg = alpha(mdof).*reshape( dc(:,:,mdof), n, ldof)./emax;
        
        %% VOLUME CONSTRAINT
        % Optionally add a volume constraint by uncommenting the following
        volfrac = 0.5;
        g = [g, sum(xPhys(:))/(n*volfrac) - 1];
        dgv = conv2(ones(nely, nelx)/(n*volfrac)./Hs,h,'same');
        dg = [dg, dgv(:)];
        
        %% DESIGN OPTIMIZATION STEP
        % Allow limited change of variables using a move limit
        sign = (x(:) - xold1).*(xold1 - xold2); % Detect oscialltions
        movelimit(sign > 0) = mlincr * movelimit(sign > 0); % Increase move limit if no oscillations
        movelimit(sign < 0) = mldecr * movelimit(sign < 0); % Decrease move limit if oscillations

        % Build approximated subproblem using reciprocal-like functions
        % subproblem = approx(x(:), f, df, g, dg);
        subproblem = approx(x(:), f, df, g, dg);
        
        xmin = max(0, x(:) - movelimit); % Set subproblem lower bound
        xmax = min(1, x(:) + movelimit); % Set subproblem upper bound
        
        % Use fmincon (interior-point algorithm) 
        % to solve the strictly convex subproblem
        options = optimoptions( 'fmincon','Display','notify-detailed', ...
        'SpecifyObjectiveGradient', true, ...
        'SpecifyConstraintGradient', true, ...
        'Algorithm','interior-point', ...
        'Hessianfcn',@subproblem.hessian);
    
        [xnew,~,~,~,lambda] = fmincon(@subproblem.objective, x(:), [], [], [], [], ...
            xmin, xmax, @subproblem.constraints, options);
        
        
        %% ORIGINAL MMA
        % For use of the original MMA (Svanberg 1987, available upon request by email) 
        % comment the "DESIGN OPTIMIZATION STEP" section and 
        % call the mmasub subroutine via
        %
        % [xnew, ~,~,~,~,~,~,~,~, low, upp] = ...
        %   mmasub(m, n, iter, x(:), zeros(n,1), ones(n,1), xold1, xold2, ...
        %   f, df, zeros(n,1), g', dg', zeros(m,n), low, upp, ...
        %   1, zeros(m,1), 1000*ones(m,1), zeros(m,1))
        %
        % Do not forget to initialize low = zeros(n,1), upp = ones(n,1)
        % outside of the optimization loop
        
        %% TERMINATION CRITERIA
        uv = (0.01 < xnew) & (xnew < 0.99); % Determination of unbounded variables 
        kktnorm = norm(df(uv) + sum(lambda.ineqnonlin'.*dg(uv, :), 2)); % Calculation of KKT norm
        deltax = mean(abs(xnew(:) - x(:))); % Average variable change
        
        %% VARIABLE UPDATE
        xold2 = xold1;
        xold1 = x(:);
        x = xnew;
        
        %% PRINT RESULTS theth
        fprintf(' %3i  f: %2.2e ', iter, -f); % Print objective
        fprintf(' g: %+2.2e ', g); % Print constraints
        fprintf(' dx: %2.2e  kkt: %2.2e ', deltax, kktnorm);
        fprintf(' E: [%5.2e, %5.2e, %5.2e] \n', strain_energy(1), strain_energy(2), strain_energy(3)); % Print strain energies
        
        %% PLOT DESIGN
        colormap(gray); imagesc(1-xPhys); caxis([0 1]); % Print design
        axis equal; axis off; drawnow;
        
        %% TERMINATION
        if (deltax < toldx) && (kktnorm < tolkkt) && prod(g < (0 + tolg))% Termination criteria
            fprintf('\n Optimization problem successfully solved \n');
            break
        end
    end
end

% Copyright (c) 2021 STIJN KOPPEN
%
% This Matlab code is written by:
% S. Koppen, M. Langelaar and F. van Keulen (June 2021)
% Delft University of Technology, Delft, Mekelweg 2, 2628 CD, 
% The Netherlands
% Faculty of Mechanical, Maritime and Materials Engineering (3ME)
% Department of Precision and Microsystems Engineering (PME)
% Specialisation in Structural Optimization and Mechanics (SOM)
% 
% Please sent your comments to: s.koppen@tudelft.nl
% or ask for a pull request on github.com/artofscience/flexures
%
% This code is intended for educational purposes and theoretical details
% are discussed in the paper
% "A simple and versatile topology optimization formulation for 
% flexure design" (2021) S. Koppen, M. Langelaar, F. Van Keulen.
% Struct Multidisc Optim
%
% Disclaimer:
% The authors reserves all rights but do not guaranty that the code is
% free from errors. Furthermore, we shall not be liable in any event
% caused by the use of the program.   
% 
% This code is based upon top88.m. The code as well as 
% a postscript version of the paper can be downloaded from the web-site: 
% http://www.topopt.dtu.dk
% "Efficient topology optimization in MATLAB using 88 lines of code (2011)"
% E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund
% Struct Multidisc Optim, Volume 43, Issue 1, p.1 - 16
