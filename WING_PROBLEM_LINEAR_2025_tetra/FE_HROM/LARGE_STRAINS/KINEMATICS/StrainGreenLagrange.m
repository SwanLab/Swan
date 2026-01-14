function Ev = StrainGreenLagrange(Fv,ndim)
% ----------------------------------------------------------------------------------------
% FUNCTION: StrainGreenLagrange
% ----------------------------------------------------------------------------------------
% PURPOSE:
%   Computes the Green–Lagrange strain tensor from the deformation gradient F, for 2D or 3D
%   configurations, using Voigt notation. This strain measure is suitable for finite strain
%   formulations, capturing geometric nonlinearities in a Lagrangian frame.
%
%   The Green–Lagrange strain is defined as:
%       E = 0.5 * (FᵗF - I)
%
%   In Voigt notation:
%     - 2D case: E = [E11; E22; 2*E12]
%     - 3D case: E = [E11; E22; E33; 2*E12; 2*E13; 2*E23]
%
% USAGE:
%   Ev = StrainGreenLagrange(Fv, ndim)
%
% INPUTS:
%   - Fv    : [ndofF × nsnap] matrix of deformation gradients at Gauss points in stacked Voigt format
%   - ndim  : Spatial dimension (2 or 3)
%
% OUTPUT:
%   - Ev    : [nstrain × nsnap] matrix of Green–Lagrange strain vectors in Voigt notation
%
% NOTES:
%   - The function is vectorized over all integration points and time steps (or snapshots).
%   - For 2D:
%       F = [F11 F12
%            F21 F22]
%       is mapped into Fv = [F11; F12; F21; F22] stacked for all Gauss points
%   - For 3D:
%       F = [F11 F12 F13
%            F21 F22 F23
%            F31 F32 F33]
%       is mapped into Fv = [F11; F12; ...; F33]
%
% REFERENCES:
%   - See /LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/README_RigidBodyMotions.pdf, page 14
%   - Symbolic derivation in /FE_HROM/LARGE_STRAINS/SYMBOLIC/FindEgreenl.m
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Date: May 2024
%   Comments by ChatGPT4, 12-May-2025
%
% SEE ALSO:
%   - StrainGreenLagrange_small
%   - PK2stress_Constitutive_Model
%   - StressesFromDisplacementsVARincre
%
% ----------------------------------------------------------------------------------------

% Green-Lagrange strain tensor 
% See
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 14
% See also /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/FindEgreenl.m

if ndim == 2
    % E = 0.5*(F'*F - I )
    % --------------------------------------------------------
% GREEN–LAGRANGE STRAIN COMPUTATION IN 2D (Voigt format)
% --------------------------------------------------------

% Number of components of the deformation gradient F in 2D (F11, F12, F21, F22)
    nF = 4; 
    % Total number of Gauss points (each with 4 F components stacked vertically)
    nelem_ngaus = size(Fv,1)/nF  ; 
    
% Index arrays to extract specific components of Fv
% Fv is assumed to be a (4×ngaus) × nsnap matrix, where rows are:
% [F11; F12; F21; F22] for each Gauss point
    FCOLS_1 = 1:4:size(Fv,1) ;
    FCOLS_2 = 2:4:size(Fv,1) ;
    FCOLS_3 = 3:4:size(Fv,1) ;
    FCOLS_4 = 4:4:size(Fv,1) ;
    %
    % Number of independent strain components in 2D Voigt notation: [E11; E22; 2E12]
    nstrain = 3; 
    
% Allocate memory for the strain tensor in Voigt format:
% one column per snapshot or time step
    Ev = zeros(nstrain*nelem_ngaus,size(Fv,2))  ;
    % Index arrays to assign values to Voigt components within Ev
    ECOLS_1 = 1:nstrain:size(Ev,1) ;
    ECOLS_2 = 2:nstrain:size(Ev,1) ;
    ECOLS_3 = 3:nstrain:size(Ev,1) ;    
    
    % E1 = E(1,1) = F1^2/2 + F4^2/2 - 1/2 ;
    Ev(ECOLS_1,:) = 0.5*(Fv(FCOLS_1,:).^2 + Fv(FCOLS_4,:).^2 -1) ;
    %  E2  = E(2,2 )=  F2^2/2 + F3^2/2 - 1/2
    Ev(ECOLS_2,:) = 0.5*(Fv(FCOLS_2,:).^2 + Fv(FCOLS_3,:).^2 -1) ;
    % E3  = E(1,2) = E(2,1) = 2* ((F1*F3)/2 + (F2*F4)/2)  % Voigt notation, no
    Ev(ECOLS_3,:) =  Fv(FCOLS_1,:).* Fv(FCOLS_3,:) +  Fv(FCOLS_2,:).* Fv(FCOLS_4,:)   ;
    
else
     
    nF = 9 ;
      nelem_ngaus = size(Fv,1)/nF  ; 
    FCOLS = cell(1,nF) ;
    for icols  =1:nF
        FCOLS{icols} = icols:nF:size(Fv,1) ;
    end
    nstrain = 6 ;
    Ev = zeros(nstrain*nelem_ngaus,size(Fv,2))  ;
    ECOLS = cell(1,nstrain) ;
    for icols = 1:nstrain
        ECOLS{icols} = icols:nstrain:size(Ev,1) ;
    end
    
    % E1 = F1^2/2 + F8^2/2 + F9^2/2 - 1/2
    Ev(ECOLS{1},:) = 0.5*(Fv(FCOLS{1},:).^2 + Fv(FCOLS{8},:).^2 + Fv(FCOLS{9},:).^2 -1) ;
    % E2 = F2^2/2 + F6^2/2 + F7^2/2 - 1/2
    Ev(ECOLS{2},:) = 0.5*(Fv(FCOLS{2},:).^2 + Fv(FCOLS{6},:).^2 + Fv(FCOLS{7},:).^2 -1) ;
    % E3 =  F3^2/2 + F4^2/2 + F5^2/2 - 1/2
    Ev(ECOLS{3},:) = 0.5*(Fv(FCOLS{3},:).^2 + Fv(FCOLS{4},:).^2 + Fv(FCOLS{5},:).^2 -1) ;
    % E4 = F2*F4 + F3*F7 + F5*F6
    Ev(ECOLS{4},:) =  Fv(FCOLS{2},:).*Fv(FCOLS{4},:) + Fv(FCOLS{3},:).*Fv(FCOLS{7},:) + Fv(FCOLS{5},:).*Fv(FCOLS{6},:)  ; 
    % E5 = F1*F5 + F3*F8 + F4*F9
    Ev(ECOLS{5},:) =  Fv(FCOLS{1},:).*Fv(FCOLS{5},:) + Fv(FCOLS{3},:).*Fv(FCOLS{8},:) + Fv(FCOLS{4},:).*Fv(FCOLS{9},:)  ; 
    % E6 = F1*F6 + F2*F9 + F7*F8
    Ev(ECOLS{6},:) =  Fv(FCOLS{1},:).*Fv(FCOLS{6},:) + Fv(FCOLS{2},:).*Fv(FCOLS{9},:) + Fv(FCOLS{7},:).*Fv(FCOLS{8},:)  ;
end