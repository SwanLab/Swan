function VAR = KineticAndStrainEnergyGet_COROT(VAR,OPERFE,DATA,MATPRO,FgradST,dCqLOC,D_QrotALL)
% Adaptation of KineticAndStrainEnergyGet.m to the corotational approach.
% % JAHO, 29-OCT-2024, HONEST GREENS, PEDRALBES CENTER, BARCELONA
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
if nargin == 0
    load('tmp.mat')
end

% Kinetic energy

if DATA.ISDYNAMIC == 1
    VAR.KINETIC_ENERGY = 0.5*(VAR.VEL'*(OPERFE.M*VAR.VEL)) ;
end

% Strain energy

VAR.STRAIN_ENERGY = 0 ;

%if isempty(OPERFE.KinternalFORCES_given)
% STANDARD OPTIONS
% ----------------
if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
    error('Option not implemented')
    for istrain = 1:DATA.MESH.nstrain
        ROWS = istrain:DATA.MESH.nstrain:length(VAR.PK2STRESS) ;
        %    VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs) ;
        VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs(:,1)) ; % Amendment 27-Jan-2024
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/02_ELASTODYNAMICS.mlx
    end
    VAR.STRAIN_ENERGY = 0.5*VAR.STRAIN_ENERGY ;
    % ----------------------------------------------------------------------------------------
else
    % JAHO, 4-May-2024
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.mlx
    % Internal forces are calculated as (see
    %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/KINEMATICS/InternalForces.m
    % FINT = Kelas*d + B^T W stress_inelastic
    
    
    % \begin{equation}
    %  \mathcal{E}^{(e)} = \dfrac{1}{2} \Par{  (\Delta \dCqLOC)^T   \DiagC{\KcLINloc}   \Delta \dCqLOC  + (\FgradFlocST-\identST)^T \DiagC{\Wecm} (\PoneFlocST - \PoneFlocLINst)      }
    % \end{equation}
    %
    % where
    %
    % \begin{equation}
    %  \Delta \dCqLOC =  \DiagC{\QrotALL}^T  \LboolCall \dC -  \dCqLOC
    % \end{equation}
    
    % FgradST - IDENTITY (gradient of displacements)
    
    for istrain = 1:DATA.MESH.ndim
        ROWS = istrain:DATA.MESH.ndim^2:length(VAR.PK1STRESS_incre) ;
        FgradST(ROWS) =FgradST(ROWS)-1 ;
    end
    
    
    VAR.Delta_dCqLOC = D_QrotALL'*(OPERFE.LboolCall*VAR.DISP) - dCqLOC ;
    VAR.STRAIN_ENERGY  = 0.5*VAR.Delta_dCqLOC'*(OPERFE.D_KcLINloc*VAR.Delta_dCqLOC);
    VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY +  FgradST'*(OPERFE.D_Wecm*VAR.PK1STRESS_incre) ;
    
 
    
end
%
% else
%
%
%     % This was implemented when testing the Dynamic Mode
%     % Decomposition, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DynamicModeDecomposition/01_ForcedVibration_2D/RunDMD_aux.mlx
%
%     if ~isstruct(OPERFE.KinternalFORCES_given)
%         K = OPERFE.KinternalFORCES_given ;
%     else
%         iMATRIX = OPERFE.KinternalFORCES_given.INTERVALS(DATA.istep) ;
%         K = OPERFE.KinternalFORCES_given.MATRICES{iMATRIX} ;
%     end
%
%     VAR.STRAIN_ENERGY = 0.5*VAR.DISP'*(K*VAR.DISP) ;
%
%
% end



% Potential energy
% -------------------------------
if isfield(DATA,'InitialPotentialEnergy')
    VAR.POTENTIAL_ENERGY = sum(OPERFE.Nst_potential_energy*VAR.DISP)  ;
    VAR.POTENTIAL_ENERGY = VAR.POTENTIAL_ENERGY + DATA.InitialPotentialEnergy ;
end

%  idim = DATA.idim_gravity; %
%  u = OPERFE.Nst*VAR.DISP;   % Displacement at all gauss points
%  %P = abs(DATA.vGRAVITY(idim))*sum(u(idim:DATA.MESH.ndim:end).*)
%
%  %VAR.POTENTIAL_ENERGY =