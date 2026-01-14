function VAR = KineticAndStrainEnergyGet(VAR,OPERFE,DATA,FgradST)

if nargin == 0
    load('tmp.mat')
end

% Kinetic energy

if DATA.ISDYNAMIC == 1
    VAR.KINETIC_ENERGY = 0.5*(VAR.VEL'*(OPERFE.M*VAR.VEL)) ;
end

% Strain energy

VAR.STRAIN_ENERGY = 0 ;

if isempty(OPERFE.KinternalFORCES_given)
    % STANDARD OPTIONS
    % ----------------
    if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0        
       % THIS IS ALSO FUNDAMENTALLY WRONG. IT IS ONLY USEFUL FOR ELASTIC
       % PROBLEMS IN WHICH ROTATIONS MAY BE LARGE, BUT STRAINS ARE SMALL  
       % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/11_AUXETIC_2Dcomp.mlx
        for istrain = 1:DATA.MESH.nstrain
            ROWS = istrain:DATA.MESH.nstrain:length(VAR.PK2STRESS) ;
            %    VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs) ;
            VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs(:,1)) ; % Amendment 27-Jan-2024
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/02_ELASTODYNAMICS.mlx
        end
        VAR.STRAIN_ENERGY = 0.5*VAR.STRAIN_ENERGY ;
        % ----------------------------------------------------------------------------------------
    else
        
        
        if DATA.SMALL_STRAIN_KINEMATICS == 1
            
            % JAHO, 4-May-2024
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.mlx
            % Internal forces are calculated as (see
            %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/KINEMATICS/InternalForces.m
            % FINT = Kelas*d + B^T W stress_inelastic
            
            for istrain = 1:DATA.MESH.nstrain
                ROWS = istrain:DATA.MESH.nstrain:length(VAR.PK2STRESS) ;
                %    VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs) ;
                VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS_incre(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs(:,1)) ;
            end
            VAR.STRAIN_ENERGY = 0.5*(VAR.STRAIN_ENERGY + VAR.DISP'*(OPERFE.KstiffLINEAR*VAR.DISP)) ;
            
        else
            % Change introduced 9-Nov-2024, KW:ST
            % see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/03_STRAIN_ROT_DIRICH.mlx
            % OBSERVATION 30-APRIL-2025: tHIS WAY OF COMPUTING STRAIN
            % ENERGY IS FUNDAMENTALLY WRONG!!! THIS IS BY NO MEANS TRUE: sum((VAR.PK1STRESS_incre(ROWS).*(FgradST(ROWS)-1)).*OPERFE.wSTs(:,1)) 
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/11_AUXETIC_2Dcomp.mlx
            % It only delivers accurate estimate in cases in which strains
            % are small/moderate
            for istrain = 1:DATA.MESH.ndim
                ROWS = istrain:DATA.MESH.ndim^2:length(VAR.PK1STRESS_incre) ;
                %    VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs) ;
                VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK1STRESS_incre(ROWS).*(FgradST(ROWS)-1)).*OPERFE.wSTs(:,1)) ;
            end
            
            for istrain = (DATA.MESH.ndim+1):DATA.MESH.ndim^2
                ROWS = istrain:DATA.MESH.ndim^2:length(VAR.PK1STRESS_incre) ;
                %    VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK2STRESS(ROWS).*VAR.GLSTRAINS(ROWS)).*OPERFE.wSTs) ;
                VAR.STRAIN_ENERGY = VAR.STRAIN_ENERGY  + sum((VAR.PK1STRESS_incre(ROWS).*FgradST(ROWS)).*OPERFE.wSTs(:,1)) ;
            end
            VAR.STRAIN_ENERGY = 0.5*(VAR.STRAIN_ENERGY + VAR.DISP'*(OPERFE.KstiffLINEAR*VAR.DISP)) ;
        end
        
        
    end
    
else
    
    
    % This was implemented when testing the Dynamic Mode
    % Decomposition, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DynamicModeDecomposition/01_ForcedVibration_2D/RunDMD_aux.mlx
    
    if ~isstruct(OPERFE.KinternalFORCES_given)
        K = OPERFE.KinternalFORCES_given ;
    else
        iMATRIX = OPERFE.KinternalFORCES_given.INTERVALS(DATA.istep) ;
        K = OPERFE.KinternalFORCES_given.MATRICES{iMATRIX} ;
    end
    
    VAR.STRAIN_ENERGY = 0.5*VAR.DISP'*(K*VAR.DISP) ;
    
    
end



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