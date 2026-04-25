% Copyright (c) 2021 STIJN KOPPEN
%
% class approx
%
% Author1:  Stijn Koppen (s.koppen@tudelft.nl)
% Author2:  Matthijs Langelaar
% Author3:  Fred van Keulen
% Date:     1 June 2021
%
% Class:        approx
% Description:  Object that builds a local approximation using
%               reciprocal-like functions
% Properties:   n   int             Number of design variables
%               m   int             Number of constraints
%               x0  float   (n x 1) Current design
%               f   float   (1 x 1) Objective value at x0
%               g   float   (1 x m) Constraint(s) values at x0
%               df  float   (n x 1) Objective sensitivities at x0
%               dg  float   (n x m) Constraint(s) sensitivities at x0
%               
%
% Examples of usage:
%
% Given a design x0 with:
% f: objective at x0
% g: constraints at x0
% df: first order objective sensitivities at x0
% dg: first order constraint sensitivities at x0
% 
% Initialse an approximated problem via, e.g.
% subproblem = approx(x0, f, df, g, dg)
% 
% Calculate approximated function values and sensitivities via:
% [fapprox, dfapprox] = subproblem.objective(x) (objective)
% [gapprox, ~, dgapprox, ~] = subproblem.constraints(x) (constraints)
%
% Diagonal Hessian of the Lagrangian is obtained via:
% h = hessian(x, lambda)


% f(y) = g(y0) + sum_i( dg/dy_i * (y_i - y0_i) ) = 
% [ g(y0) - sum_i( dg/dy_i * y0_i ) ] + sum_i( dg/dy_i * y_i )
%
% f(x) = g(x0) + sum_i( dg/dx_i * dx_i/dy_i * (y_i(x_i) - y0_i(x0_i)) ) =
% [g(x0) - sum_i( dg/dx_i * dx_i/dy_i * y0_i(x0_i) ) ] + 
% sum_i( dg/dx_i * dx_i/dy_i * y_i(x_i) )

% df(x)/dx_i = dg/dy_i * dy_i/dx_i(x)
% d2f(x)/d2x_i = dg/dy_i * d2y_i/d2x_i(x)

% for linear:
% y = x
% dy/dx = 1
% d2y/d2x = 0

% x = y
% dx/dy = 1
% d2x/d2y = 0

% f(x) = g(x0) + sum_i( dg/dx_i * dx_i/dy_i * (y_i(x_i) - y0_i(x0_i)) )
% f(x) = g(x0) + sum_i( dg/dx_i * (x_i - x0_i) )

% for mma:  
% y =           (u - x)^-1
% dy/dx =       (u - x)^-2
% d2y/d2x =     2*(u - x)^-3

% dg/dy = dg/dx (provided) * dx/dy
% x =           u -y^-1
% dx/dy =       y^-2 = (u - x)^2
% d2x/d2y =     -2*y^-3 = -2*(u - x)^3


classdef approx
    properties
        x0; n; f; g; df; dg; % Inputs (description above)
        f0; % Approximated objective value at x = 0
        g0; % Approximated constraint values at x = 0
        low = -1; % Lower asymptote, note low < xmin
        upp = 2; % Upper asymptote, note upp > xmax
        dfdy; % Derivative of objective wrt intervening variable y = 1/(x - low)
        dgdy; % Derivative of constraints wrt intervening variables y = 1/(upp - x)
        dx0l; % x0 - low
        dx0u; % upp - x0
    end
    methods
        function obj = approx(x0, f, df, g, dg)
            obj.x0 = x0; 
            obj.f = f; obj.g = g; 
            obj.df = df; obj.dg = dg; % Assign inputs
            
            obj.n = length(x0);
            
            obj.dx0l = obj.x0  - obj.low;
            obj.dx0u = obj.upp - obj.x0;
            
            
            obj.dfdy = -obj.df.*obj.dx0l.^2;
            obj.dgdy =  obj.dg.*obj.dx0u.^2;
            
            obj.f0 = obj.f - 1./obj.dx0l'*obj.dfdy;
            obj.g0 = obj.g - 1./obj.dx0u'*obj.dgdy;
        end
        
        function [fapprox, dfapprox] = objective(obj, x)
            dxl = x - obj.low; % Local variable dxl
            fapprox = obj.f0 + 1./dxl'*obj.dfdy;
            if nargout > 1
                dfapprox = -obj.dfdy.*1./dxl.^2;
            end
        end
        
        function [gapprox, happrox, dgapprox, dhapprox] = constraints(obj, x)
            dxu = obj.upp - x; % Local variable dxu
            gapprox = obj.g0 + (1./dxu)'*obj.dgdy; % Nonlinear inequality constraints
            happrox = []; % Nonlinear equality constraints (none)
            if nargout > 2
                dgapprox = obj.dgdy./dxu.^2; % Sensitivities of nonlinear inequality constraints
                dhapprox = []; % Sensitivities of nonlinear equality constraints (none)
            end
        end
        
        function h = hessian(obj, x, lambda)
            hf = 2*obj.dfdy./(x - obj.low).^3; % Diagonal hessian of approximated objective
            hg = 2*obj.dgdy./(obj.upp - x).^3; % Diagonal hessian of approximated constraints
            
            h = spdiags(hf, 0, obj.n, obj.n);
            
            % Hessian of Lagrangian (L = fapprox + sum_i(lambda_i * gapprox_i))
            for i=1:length(lambda.ineqnonlin) % Loop over number of constraints
                h = h + lambda.ineqnonlin(i)*spdiags(hg(:,i), 0, obj.n, obj.n);
            end       
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
