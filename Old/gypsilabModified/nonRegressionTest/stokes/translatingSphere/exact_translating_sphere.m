function [U_exact] = exact_translating_sphere(sphere_translational_velocity,sphere_radius,mesh_eval)
%+========================================================================+
%|                                                                        |
%|           This function uses the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : exact_translating_sphere.m                    |
%|    #    |   VERSION    : 0.10                                          |
%|   _#_   |   AUTHOR(S)  : Luca Berti                                    |
%|  ( # )  |   CREATION   : 25.12.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+
% COMPUTE THE EXACT VELOCITY FIELD DETERMINED BY A TRANSLATING SPHERE IN 
% AN INFINITE FLUID ALONG THE Z AXIS 
% Reference: Happel-Brenner, "Low Reynolds number hydrodynamics"
% (pag. 119 and formulas 4-17.17 and 4-17.18) transformed into Cartesian 
% coordinates (p.507 A-15.29 for the transformation)

% Input -> mesh_eval: mesh containing the evaluation points for the exact
%                     solution 
%          sphere_translational_velocity -> translational velocity of the
%                                           vertically moving sphere
%                                           [0,0,W]
% Output -> U_exact:  3xN_eval vector containing the values of the exact
%                     solution at the evaluation point, expressed IN CARTESIAN
%                     COORDINATES

W = sphere_translational_velocity(3);
% Exact solutions taken from the reference - SPHERICAL COORDINATES
u_r = @(x,y,z) -W.*cos(acos(z./sqrt(x.^2+y.^2+z.^2))).*(sphere_radius^3./(2*sqrt(x.^2+y.^2+z.^2).^3) -3*sphere_radius./(2*sqrt(x.^2+y.^2+z.^2)));
u_theta = @(x,y,z) -W.*sin(acos(z./sqrt(x.^2+y.^2+z.^2))).*(sphere_radius^3./(4*sqrt(x.^2+y.^2+z.^2).^3)+3*sphere_radius./(4*sqrt(x.^2+y.^2+z.^2)));
u_phi =@(x,y,z) 0;

X = mesh_eval.vtx;
N_eval = size(X,1);

% Evaluation vector - SPHERICAL COORDINATES
U_exact = zeros(N_eval*3,1);
U_exact(1:N_eval,1) = feval(u_r,X(:,1),X(:,2),X(:,3));
U_exact(N_eval+1:2*N_eval,1) = feval(u_theta,X(:,1),X(:,2),X(:,3));
U_exact(2*N_eval+1:end,1) = feval(u_phi,X(:,1),X(:,2),X(:,3));

% Expression of the evaluation points in spherical coordinates, to extract
% (r,theta,phi) needed to transform the exact velocities from spherical to
% Cartesian coordinates
Y = X;
[X(:,3),X(:,2),X(:,1)] = cart2sph(Y(:,1),Y(:,2),Y(:,3));
% Latitude is needed instead of colatitude (as computed by Matlab)
X(:,2) = pi/2-X(:,2);

temporary = zeros(N_eval*3,1);
% Transformation of the exact velocity vector field from spherical to
% Cartesian coordinates - CARTESIAN COORDINATES
temporary(1:N_eval,1) = U_exact(1:N_eval,1).*sin(X(:,2)).*cos(X(:,3)) + ...
                          U_exact(N_eval+1:N_eval*2,1).*cos(X(:,2)).*cos(X(:,3))- ...
                          U_exact(1+2*N_eval:end,1).*sin(X(:,3));

temporary(1+N_eval:2*N_eval,1) = U_exact(1:N_eval,1).*sin(X(:,2)).*sin(X(:,3)) + ...
                          U_exact(N_eval+1:N_eval*2,1).*cos(X(:,2)).*sin(X(:,3)) + ...
                          U_exact(1+2*N_eval:end,1).*cos(X(:,3));
                       
temporary(1+2*N_eval:end,1) = U_exact(1:N_eval,1).*cos(X(:,2)) - ...
                          U_exact(N_eval+1:N_eval*2,1).*sin(X(:,2));
U_exact = temporary; 


end