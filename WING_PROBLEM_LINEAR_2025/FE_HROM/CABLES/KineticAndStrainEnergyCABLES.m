function VAR = KineticAndStrainEnergyCABLES(VAR,OPERFE,DATA,MATPRO)

if nargin == 0
    load('tmp.mat')
end

% Kinetic energy


VAR.KINETIC_ENERGY = 0.5*(VAR.VEL'*(OPERFE.M*VAR.VEL)) ;

% Strain energy

%VAR.STRAIN_ENERGY = 0 ;

%if isempty(OPERFE.KinternalFORCES_given)
    %for istrain = 1:DATA.MESH.nstrain
     %   ROWS = istrain:DATA.MESH.nstrain:length(VAR.PK2STRESS) ;
        VAR.STRAIN_ENERGY =  0.5*sum((VAR.TENSION.*VAR.STRAIN).*OPERFE.wSTs);
    %end
  %  VAR.STRAIN_ENERGY = 0.5*VAR.STRAIN_ENERGY ;
    
%else
    
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
%     
%     
%     
%     VAR.STRAIN_ENERGY = 0.5*VAR.DISP'*(K*VAR.DISP) ;
%     
%     
%     
%     
%     
% end
% 
% 
% 
% % Potential energy
% % -------------------------------
% VAR.POTENTIAL_ENERGY = sum(OPERFE.Nst_potential_energy*VAR.DISP)  ;
% VAR.POTENTIAL_ENERGY = VAR.POTENTIAL_ENERGY + DATA.InitialPotentialEnergy ;
% 
% %  idim = DATA.idim_gravity; %
% %  u = OPERFE.Nst*VAR.DISP;   % Displacement at all gauss points
% %  %P = abs(DATA.vGRAVITY(idim))*sum(u(idim:DATA.MESH.ndim:end).*)
% %
% %  %VAR.POTENTIAL_ENERGY =