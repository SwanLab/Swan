function VAR =  UpdateVelocityAcceleration(VAR,TIMEPARAM,DeltaT)
%--------------------------------------------------------------------------
% VAR = UpdateVelocityAcceleration(VAR, TIMEPARAM, DeltaT)
%
% PURPOSE:
%   Updates the velocity and acceleration fields in VAR using the
%   Newmark–Beta–Gamma time integration scheme (specifically, the 
%   Newmark-Bossak formulation) after solving for displacements.
%
% DESCRIPTION:
%   This function computes the current velocity and acceleration
%   vectors at time step n+1 based on the known displacements and
%   predictor values from the previous time step.
%
%   The formulas implemented are:
%
%     v^{n+1} = (γ / (β Δt)) * d^{n+1} + [vTILDE - (γ / (β Δt)) * dTILDE]
%     a^{n+1} = (1 / (β Δt²)) * (d^{n+1} - dTILDE)
%
%   These are derived from the Newmark-Bossak method, which introduces
%   numerical damping to stabilize time integration in structural dynamics.
%
% INPUT:
%   VAR       : Structure containing current and predictor values of
%               displacement, velocity, and acceleration.
%               Fields used: DISP, VELtilde, DISPtilde
%
%   TIMEPARAM : Structure containing Newmark–Bossak parameters:
%               - TYPE  : must be 1 (for Bossak scheme)
%               - BETA  : β parameter (displacement weighting)
%               - GAMMA : γ parameter (velocity weighting)
%
%   DeltaT    : Time step size Δt
%
% OUTPUT:
%   VAR       : Updated structure, including new velocity (VAR.VEL)
%               and acceleration (VAR.ACEL)
%
% REFERENCES:
%   - Newmark, N. M. (1959). A method of computation for structural dynamics.
%   - Bossak, M. (1983). A method for the direct integration of structural dynamics.
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp2.mat')
end






if TIMEPARAM.TYPE==1
    % Newmark-Bossak
    
    %*************
    % Velocity
    %************
    
    % %  \begin{equation*}
    % \label{eq:vel3}
    % \begin{split}
    %  \v^{n+1} & =    \dfrac{\gamma}{\beta \Delta t_{n+1}}  \d^{n+1}    + \Par{\vTIL^{n+1} - \dfrac{\gamma}{\beta \Delta t_{n+1}}  \dTIL^{n+1}}
    %  \end{split}
    % \end{equation*}
    
    
    VAR.VEL =  TIMEPARAM.GAMMA/(TIMEPARAM.BETA*DeltaT)*VAR.DISP + (VAR.VELtilde - TIMEPARAM.GAMMA/(TIMEPARAM.BETA*DeltaT)*VAR.DISPtilde) ;
    
    
    %*************
    % ACCELERATION
    %************
    %     \begin{equation*}
    % \label{eq:42lkk}
    % \begin{split}
    %  \a^{n+1} & =  \dfrac{1}{\beta \Delta t_{n+1}^2} \Par{\d^{n+1} -  \dTIL^{n+1}}
    %  \end{split}
    % \end{equation*}
    
    VAR.ACEL = 1/(TIMEPARAM.BETA*DeltaT^2)*(VAR.DISP-VAR.DISPtilde) ; 
    
    
end