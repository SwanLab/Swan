function [DATA,VAR_n,SNAP] = INITIALIZATIONvar_BUBFE(DATA,INICOND)
if nargin ==1
    INICOND = [] ;
end
% Initialization displacement, velocity and acceleration, on the one hand
% and other variables defined at Gauss points
% JAHO 26-Nov-2020
if nargin == 0
    load('tmp3.mat')
end


if isfield(DATA.STORE,'TOLERANCE_SVD_COMPRESSION') &&  isnumeric(DATA.STORE.TOLERANCE_SVD_COMPRESSION)
    fff = fieldnames(DATA.STORE.VAR) ;
    TOL = DATA.STORE.TOLERANCE_SVD_COMPRESSION ;
    DATA.STORE.TOLERANCE_SVD_COMPRESSION = [];
    for  iii = 1:length(fff)
        DATA.STORE.TOLERANCE_SVD_COMPRESSION.(fff{iii}) = TOL ;
    end
end


DATA = DefaultField(DATA,'INTEGRATION_SCHEME',[]) ;
DATA.INTEGRATION_SCHEME = DefaultField(DATA.INTEGRATION_SCHEME,'steps_KINEM_VARIABLES_TO_STORE',1) ;
DATA.INTEGRATION_SCHEME = DefaultField(DATA.INTEGRATION_SCHEME,'steps_STRESS_VARIABLES_TO_STORE',1) ;
% Number of previous time steps information
% to be stored  (with a view towards extending to other time-integrating systems
% (or for doing extrapolations using, say, the SVD...). By default = 1 )

nstoreKINEM = DATA.INTEGRATION_SCHEME.steps_KINEM_VARIABLES_TO_STORE  ;
nstoreSTRESS = DATA.INTEGRATION_SCHEME.steps_STRESS_VARIABLES_TO_STORE  ;
% --------------------
% Kinematic variables
% --------------------
ndof = DATA.MESH.ndof ;

if isempty(INICOND) || ~isfield(INICOND,'DISP')
    INICOND.DISP = zeros(ndof,1) ;
    
end

VAR_n.DISP = repmat(INICOND.DISP,1,nstoreKINEM) ;  % displacements, initial conditions in
VAR_n.BUB_DISP = VAR_n.DISP ; 
VAR_n.ELAST_DISP = VAR_n.DISP ; 

VAR_n.RESID = zeros(ndof,1)  ;
VAR_n.FEXT = zeros(ndof,1)  ;
VAR_n.FINT = zeros(ndof,1)  ;



DATA = DefaultField(DATA,'ISDYNAMIC',0) ;
VAR_n.STRAIN_ENERGY = 0 ;
if DATA.ISDYNAMIC==1
    VAR_n.VEL  =repmat(INICOND.VELOC,1,nstoreKINEM)  ; % Velocity, initial conditions in
    VAR_n.ACEL = zeros(ndof,nstoreKINEM)  ;             % Acceleration
    
    VAR_n.KINETIC_ENERGY = 0 ;
   % VAR_n.STRAIN_ENERGY = 0 ;
    
    % Initial
    DATA = DefaultField(DATA,'InitialPotentialEnergy',0) ; 
    VAR_n.POTENTIAL_ENERGY = DATA.InitialPotentialEnergy ;
end

DATA = DefaultField(DATA,'FOLLOWER_LOADS',[]) ;
if ~isempty(DATA.FOLLOWER_LOADS)
    DATA.FOLLOWER_LOADS = DefaultField(DATA.FOLLOWER_LOADS,'HYDROSTATIC',[]) ;
    if ~isempty(DATA.FOLLOWER_LOADS.HYDROSTATIC)
        VAR_n.FEXTpress= zeros(ndof,1)  ;
    end
end




% --------------------------
% Gauss variables
%------------------------
VAR_n.PK2STRESS = zeros(DATA.MESH.ndofSTRESS,1,nstoreSTRESS) ;

ndofPK1stress = DATA.MESH.ngausT*DATA.MESH.ndim^2 ;

%if DATA.SMALL_STRAIN_KINEMATICS ==1 && DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1
%else
    VAR_n.PK1STRESS = zeros(ndofPK1stress,1,nstoreSTRESS) ;
%end

VAR_n.GLSTRAINS = zeros(DATA.MESH.ndofSTRESS,1,nstoreSTRESS) ;

VAR_n.CAUCHY_STRESS = zeros(DATA.MESH.ndofSTRESS,1,nstoreSTRESS) ;

VAR_n.VONMISES_CAUCHY_STRESS = zeros(DATA.MESH.ngausT,1,nstoreSTRESS) ;

if ~isempty(DATA.FOLLOWER_LOADS)
    if ~isempty(DATA.FOLLOWER_LOADS.HYDROSTATIC)
        VAR_n.tPRESS= zeros(DATA.MESH.HYDRO.ngausT*DATA.MESH.ndim,1)  ;
    end
end

% J2 model variables
% ------------------
switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
    case 'SMALL_STRAINS_J2_PLASTICITY'
        VAR_n.PlasticStrains =  repmat(INICOND.PlasticStrains,1,nstoreSTRESS) ; ;  % Plastic strains
        VAR_n.InternalVarStrain = repmat(INICOND.InternalVarStrain,1,nstoreSTRESS) ;  % STRAIN-LIKE INTERNAL VARIABLE strains
        VAR_n.YieldStress =  repmat(INICOND.YieldStress,1,nstoreSTRESS) ;  % YIELD STRESS/STRESS-LIKE INTERNAL VARIABLE strains
end
 


DATA = DefaultField(DATA,'SKIP_PART_STORE',0) ;
DATA = DefaultField(DATA,'ENABLE_PRINTING_NONCONVERGED_ITERATIONS',0) ;

SNAP = [] ;
SNAP_ITER = [] ;

if DATA.SKIP_PART_STORE == 0
    % INITIZALIZATION SNAPSHOTS
    % -------------------------
    fff =fieldnames(DATA.STORE.VAR) ;
    iclusterSTORE = 1;
    
    nstepsLOC = DATA.STORE.NSTEPS_CLUSTER{iclusterSTORE} ;
    for iii=1:length(fff)
        if DATA.STORE.VAR.(fff{iii})
            nvar = size(VAR_n.(fff{iii}),1) ;
            SNAP.(fff{iii}) = zeros(nvar,length(nstepsLOC))  ;
            SNAP.(fff{iii})(:,1) = VAR_n.(fff{iii})(:,end) ;
        end
    end
    if  DATA.ENABLE_PRINTING_NONCONVERGED_ITERATIONS == 1
        nstepsITER= DATA.NEWTON_RAPHSON.NMAXiter + 1;
        for iii=1:length(fff)
            if DATA.STORE.VAR.(fff{iii})
                nvar = size(VAR_n.(fff{iii}),1) ;
                SNAP_ITER.(fff{iii}) = zeros(nvar,nstepsITER)  ;
                % SNAP_ITER.(fff{iii})(:,1) = VAR_n.(fff{iii})(:,end) ;
            end
        end
        DATA.SNAP_ITER = SNAP_ITER ;
    else
        DATA.SNAP_ITER = [] ;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VARIABLES TO BE PRINTED
    %-----------------------------------------------------------------
    NODESV_PROP = DATA.PRINT.NODESV_PROP ;
    GAUSSV_PROP = DATA.PRINT.GAUSSV_PROP ;
    
    
    NODESV_PROP.DISP.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.DISP.LEGEND = 'DISPLAC.' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.DISP.COMP = {'u-x','u-y'} ;
    else
        NODESV_PROP.DISP.COMP = {'u-x','u-y','u-z'} ;
    end
    
      NODESV_PROP.BUB_DISP.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.BUB_DISP.LEGEND = 'BUBBLE DISPLAC.' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.BUB_DISP.COMP = {'uBUB-x','uBUB-y'} ;
    else
        NODESV_PROP.BUB_DISP.COMP = {'uBUB-x','uBUB-y','uBUB-z'} ;
    end
    
     NODESV_PROP.ELAST_DISP.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.ELAST_DISP.LEGEND = 'ELAST DISPLAC.' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.ELAST_DISP.COMP = {'uEL-x','uEL-y'} ;
    else
        NODESV_PROP.ELAST_DISP.COMP = {'uEL-x','uEL-y','uEL-z'} ;
    end
    
     
    
    
    
    
    
    NODESV_PROP.RESID.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.RESID.LEGEND = 'RESIDUAL' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.RESID.COMP = {'R-x','R-y'} ;
    else
        NODESV_PROP.RESID.COMP = {'R-x','R-y','R-z'} ;
    end
    
    NODESV_PROP.FEXT.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.FEXT.LEGEND = 'FEXT' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.FEXT.COMP = {'F-x','F-y'} ;
    else
        NODESV_PROP.FEXT.COMP = {'F-x','F-y','F-z'} ;
    end
    
    
    NODESV_PROP.FEXTpress.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.FEXTpress.LEGEND = 'FEXTpress' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.FEXTpress.COMP = {'F-x','F-y'} ;
    else
        NODESV_PROP.FEXTpress.COMP = {'F-x','F-y','F-z'} ;
    end
    
    
    
    NODESV_PROP.ACEL.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.ACEL.LEGEND = 'ACEL.' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.ACEL.COMP = {'a-x','a-y'} ;
    else
        NODESV_PROP.ACEL.COMP = {'a-x','a-y','a-z'} ;
    end
    
    
    
    NODESV_PROP.VEL.TYPE = 'Vector' ;  % For printing in GID
    NODESV_PROP.VEL.LEGEND = 'VELOC.' ;  % GID's printing
    if DATA.MESH.ndim == 2
        NODESV_PROP.VEL.COMP = {'v-x','v-y'} ;
    else
        NODESV_PROP.VEL.COMP = {'v-x','v-y','v-z'} ;
    end
    
    
    
    
    GAUSSV_PROP.PK2STRESS.LEGEND = 'PK2STRESS' ;
    
    
    if DATA.MESH.ndim  == 2
        GAUSSV_PROP.PK2STRESS.TYPE = 'Matrix' ;
        GAUSSV_PROP.PK2STRESS.COMP ={'S-X','S-Y','S-Z','S-XY'} ;
        switch DATA.typePROBLEM
            case 'pstress'
                GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1,2,4,3] ;
            otherwise
                if DATA.MESH.nstrain==3
                    warning('Stresses in the z-direction have been not computed...')
                end
                GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1,2,4,3] ;  % sxx, syy, sxy, szz
        end
        %GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1,2,3,4] ;
    elseif DATA.MESH.ndim  == 3
        GAUSSV_PROP.PK2STRESS.TYPE = 'Matrix' ;
        GAUSSV_PROP.PK2STRESS.COMP ={'S-X','S-Y','S-Z','S-XY','S-YZ','S-XZ'} ;
        GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1 2 3 6 4 5] ;
    end
    
    
    
    
    GAUSSV_PROP.CAUCHY_STRESS.LEGEND = 'CAUCHY_STRESS' ;
    
    
    if DATA.MESH.ndim   == 2  %&& DATA.MESH.nstrain == 4
        %  error('To be implemented')
        GAUSSV_PROP.CAUCHY_STRESS.TYPE = 'Matrix' ;
        GAUSSV_PROP.CAUCHY_STRESS.COMP ={'s-X','s-Y','s-Z','s-XY'} ;
        switch DATA.typePROBLEM
            case 'pstress'
                GAUSSV_PROP.CAUCHY_STRESS.CONVERSION_GID =[1,2,4,3] ;
            otherwise
                if DATA.MESH.nstrain==3
                    warning('Stresses in the z-direction have been not computed...')
                end
                GAUSSV_PROP.CAUCHY_STRESS.CONVERSION_GID =[1,2,4,3] ;
        end
    elseif DATA.MESH.ndim  == 3
        GAUSSV_PROP.CAUCHY_STRESS.TYPE = 'Matrix' ;
        GAUSSV_PROP.CAUCHY_STRESS.COMP ={'s-X','s-Y','s-Z','s-XY','s-YZ','s-XZ'} ;
        GAUSSV_PROP.CAUCHY_STRESS.CONVERSION_GID =[1 2 3 6 4 5] ;
    end
    
    %%%%%%%%%%%%%%%%%%%%
    
    GAUSSV_PROP.GLSTRAINS.LEGEND = 'GreenLGstrain' ;
    
    
    if DATA.MESH.ndim  == 2
        GAUSSV_PROP.GLSTRAINS.TYPE = 'Matrix' ;
        GAUSSV_PROP.GLSTRAINS.COMP ={'E-X','E-Y','E-Z','E-XY'} ;
        switch DATA.typePROBLEM
            case 'pstress'
                warning('Strains in the z-direction have been not computed...')
                GAUSSV_PROP.GLSTRAINS.CONVERSION_GID =[1,2,4,3] ;
            otherwise
                GAUSSV_PROP.GLSTRAINS.CONVERSION_GID =[1,2,4,3] ;
        end
    else
        GAUSSV_PROP.GLSTRAINS.TYPE = 'Matrix' ;
        GAUSSV_PROP.GLSTRAINS.COMP ={'E-X','E-Y','E-Z','E-XY','E-YZ','E-XZ'} ;
        GAUSSV_PROP.GLSTRAINS.CONVERSION_GID =[1 2 3 6 4 5] ;
    end
    
    % Cauchy stresses
    % ---------------------
    GAUSSV_PROP.VONMISES_CAUCHY_STRESS.LEGEND = 'CAUCHY VONMISES' ;
    GAUSSV_PROP.VONMISES_CAUCHY_STRESS.TYPE = 'Scalar' ;
    
    % J2 small strain variables
    % --------------------------
    switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
        case 'SMALL_STRAINS_J2_PLASTICITY'
          % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
            % PlasticStrains
            % *************************************************************************
            GAUSSV_PROP.PlasticStrains.LEGEND = 'PlasticStrains' ;
            if DATA.MESH.ndim  == 2
                GAUSSV_PROP.PlasticStrains.TYPE = 'Matrix' ;
                GAUSSV_PROP.PlasticStrains.COMP ={'EP-X','EP-Y','EP-Z','EP-XY'} ;
                switch DATA.typePROBLEM
                    case 'pstress'
                        warning('Strains in the z-direction have been not computed...')
                        GAUSSV_PROP.PlasticStrains.CONVERSION_GID =[1,2,4,3] ;
                    otherwise
                        GAUSSV_PROP.PlasticStrains.CONVERSION_GID =[1,2,4,3] ;
                end
            else
                GAUSSV_PROP.PlasticStrains.TYPE = 'Matrix' ;
                GAUSSV_PROP.PlasticStrains.COMP ={'E-X','E-Y','E-Z','E-XY','E-YZ','E-XZ'} ;
                GAUSSV_PROP.PlasticStrains.CONVERSION_GID =[1 2 3 6 4 5] ;
            end
            
            % InternalVarStrain 
            % -----------------------------
            GAUSSV_PROP.InternalVarStrain.LEGEND = 'InternalVarStrain' ;
            GAUSSV_PROP.InternalVarStrain.TYPE = 'Scalar' ;
            
            % InternalVarStrain 
            % -----------------------------
            GAUSSV_PROP.YieldStress.LEGEND = 'YieldStress' ;
            GAUSSV_PROP.YieldStress.TYPE = 'Scalar' ;
            
    end
    
    
    
    
    DATA.PRINT.NODESV_PROP = NODESV_PROP  ;
    DATA.PRINT.GAUSSV_PROP  =GAUSSV_PROP  ;
    
    
end


% % INFORMATION FOR PRINTING VARIABLES IN GID
% % ----------------------------------
% DATA.PRINT.DISP.STORE =  1 ;  % Store in the corresponding snapshot matrix
% NODESV_PROP.VEL.STORE = 0;
% NODESV_PROP.ACEL.STORE =  0 ;
%
% NODESV_PROP.DISP.PRINT = 1 ;  % GID's printing
% NODESV_PROP.VEL.PRINT = 0;  %
% NODESV_PROP.ACEL.PRINT = 0 ;  %
%
% % INFORMATION FOR PRINTING VARIABLES BY GID
% % -----------------------------------------
%
% NODESV_PROP.DISP.TYPE = 'Vector' ;  % For printing in GID
% NODESV_PROP.DISP.LEGEND = 'DISPLAC.' ;  % GID's printing
% if DATA.MESH.ndim == 2
%     NODESV_PROP.DISP.COMP = {'u-x','u-y'} ;
% else
%     NODESV_PROP.DISP.COMP = {'u-x','u-y','u-z'} ;
% end
%
% % VARIABLES DEFINED AT INTEGRATION POINTS
% % ---------------------------------------
% GAUSSV.PK2STRESS = zeros(DATA.MESH.ndofSTRESS,nstoreSTRESS) ;  % 2nd Piola-Kirchhoff stress
% GAUSSV.GLSTRAINS = zeros(DATA.MESH.ndofSTRESS,nstoreSTRESS) ;  % Green-Lagrange strains
%
% GAUSSV_PROP.PK2STRESS.PRINT = 1 ;
% GAUSSV_PROP.GLSTRAINS.PRINT = 0 ;
%
% GAUSSV_PROP.PK2STRESS.STORE = 1 ;
% GAUSSV_PROP.GLSTRAINS.STORE = 1 ;
%
%
%
% if DATA.MESH.ndim  == 2 && DATA.MESH.nstrain == 3
%     GAUSSV_PROP.PK2STRESS.TYPE = 'Vector' ;
%     GAUSSV_PROP.PK2STRESS.COMP ={'S-X','S-Y','S-XY'} ;
%     GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1,2,3] ;
% elseif DATA.MESH.ndim  == 2 && DATA.MESH.nstrain == 4
%     error('To be implemented')
%     GAUSSV_PROP.stressST.TYPE = 'Matrix' ;
%     GAUSSV_PROP.stressST.COMP ={'S-X','S-Y','S-Z','S-XY'} ;
%     GAUSSV_PROP.stressST.CONVERSION_GID =[1,2,3,4] ;
% elseif DATA.MESH.ndim  == 3
%     GAUSSV_PROP.PK2STRESS.TYPE = 'Matrix' ;
%     GAUSSV_PROP.PK2STRESS.COMP ={'S-X','S-Y','S-Z','S-XY','S-YZ','S-XZ'} ;
%     GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1 2 3 6 4 5] ;
% end
% GAUSSV_PROP.PK2STRESS.LEGEND = 'PK2STRESS' ;
%
%
%
%
%
%
%
%
%
%
% %
% % switch DATA.TypeImplementation
% %     case 'J2_plasticity_small_strains'
% %         GAUSSV_n.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
% %         GAUSSV_n.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
% %         GAUSSV_n.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
% %
% %         GAUSSV_PROP.alpha.TYPE = 'Scalar' ;
% %         GAUSSV_PROP.alpha.PRINT = 1 ;  % GID's printing
% %         GAUSSV_PROP.alpha.LEGEND = 'INT.VAR.' ;
% %
% %         GAUSSV_n.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
% %         GAUSSV_n.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
% %         GAUSSV_n.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
% %         GAUSSV_n.VonMises = zeros(OPERfe.ngaus,1);; % stressST-like internal variable
% %
% %
% %
% %
% %         GAUSSV_PROP.alpha.TYPE = 'Scalar' ;
% %
% %         DATA.PRINT_GID = DefaultField(DATA.PRINT_GID,'INTERNAL_VARIABLES',1) ;
% %         if DATA.PRINT_GID.INTERNAL_VARIABLES == 1
% %             GAUSSV_PROP.alpha.PRINT = 1 ;  % GID's printing
% %         else
% %             GAUSSV_PROP.alpha.PRINT = 0 ;
% %         end
% %         GAUSSV_PROP.alpha.LEGEND = 'INT.VAR.' ;
% %
% %
% %         DATA.STORE = DefaultField(DATA.STORE,'INTERNAL_VARIABLES',1) ;
% %         if DATA.STORE.INTERNAL_VARIABLES == 0
% %             GAUSSV_PROP.alpha.STORE = 0 ;
% %         end
% %
% %
% %         GAUSSV_PROP.VonMises.TYPE = 'Scalar' ;
% %         GAUSSV_PROP.VonMises.PRINT = 1 ;  % GID's printing
% %         GAUSSV_PROP.VonMises.LEGEND = 'Von Mises' ;
% %     otherwise
% %         error('Option not implemented')
% % end
% %
% %
% %
% % % -----------------------------------------------------------------------
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % GAUSSV_np1.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
% % % GAUSSV_np1.stressST = zeros(OPERfe.ngausT,1);  % stressSTes at time n+1
% % % GAUSSV_np1.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
% % % GAUSSV_np1.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
% % % % -----------------------------------------------------------------------
