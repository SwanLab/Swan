function [DATA_cl,VAR_n,SNAP] = INITIALIZATIONvar_cluster(DATA_cl,INICOND)
if nargin ==0
    load('tmp.mat');
end
% Initialization displacement, velocity and acceleration, on the one hand
% and other variables defined at Gauss points
% Clustering (local basis)
% JAHO 22-Feb-2021/
%
% DATA = DefaultField(DATA,'INTEGRATION_SCHEME',[]) ;
% DATA.INTEGRATION_SCHEME = DefaultField(DATA.INTEGRATION_SCHEME,'steps_KINEM_VARIABLES_TO_STORE',1) ;
% DATA.INTEGRATION_SCHEME = DefaultField(DATA.INTEGRATION_SCHEME,'steps_STRESS_VARIABLES_TO_STORE',1) ;
% Number of previous time steps information
% to be stored  (with a view towards extending to other time-integrating systems
% (or for doing extrapolations using, say, the SVD...). By default = 1 )

nstoreKINEM =1 ;% DATA.INTEGRATION_SCHEME.steps_KINEM_VARIABLES_TO_STORE  ;
nstoreSTRESS =1 ; % DATA.INTEGRATION_SCHEME.steps_STRESS_VARIABLES_TO_STORE  ;
% --------------------
% Kinematic variables
% --------------------
nclusters = length(DATA_cl.VAR) ;
ndof_cl = zeros(nclusters,1) ;
ndofSTRESS_cl = zeros(nclusters,1) ;
ngausT_cl = zeros(nclusters,1) ;

% Loop over clusters  (local basis)
for icluster = 1:nclusters
   DATA = DATA_cl.VAR{icluster} ;
    ndof_cl(icluster) = DATA.MESH.ndof ;
    ndofSTRESS_cl(icluster) = DATA.MESH.ndofSTRESS;
    ngausT_cl(icluster) = DATA.MESH.ngausT;
end
DATA = DATA_cl.COMMON  ;
ndof = max(ndof_cl) ;
ndofSTRESS = max(ndofSTRESS_cl) ;
ngausT = max(ngausT_cl) ;


if isempty(INICOND)
    INICOND.DISP = zeros(ndof,1) ;
end

VAR_n.DISP = repmat(INICOND.DISP,1,nstoreKINEM) ;  % displacements, initial conditions in
VAR_n.RESID = zeros(ndof,1)  ;             % Acceleration
VAR_n.FEXT = zeros(ndof,1)  ;             % Acceleration
% DATA = DefaultField(DATA,'ISDYNAMIC',0) ;
% if DATA.ISDYNAMIC==1
%     VAR_n.VELOC  =repmat(INICOND.VELOC,1,nstoreKINEM)  ; % Velocity, initial conditions in
%     VAR_n.ACEL = zeros(ndof,nstoreKINEM)  ;             % Acceleration
% end


% Gauss variables
VAR_n.PK2STRESS = zeros(ndofSTRESS,1,nstoreSTRESS) ;
VAR_n.GLSTRAINS = zeros(ndofSTRESS,1,nstoreSTRESS) ;

VAR_n.CAUCHY_STRESS = zeros(ndofSTRESS,1,nstoreSTRESS) ;

VAR_n.VONMISES_CAUCHY_STRESS = zeros(ngausT,1,nstoreSTRESS) ;


% INITIZALIZATION SNAPSHOTS
% -------------------------
fff =fieldnames(DATA.STORE.VAR) ;
SNAP = [] ;
iclusterSTORE = 1;

nstepsLOC = DATA.STORE.NSTEPS_CLUSTER{iclusterSTORE} ;
for iii=1:length(fff)
    if DATA.STORE.VAR.(fff{iii})
        nvar = size(VAR_n.(fff{iii}),1) ;
        SNAP.(fff{iii}) = zeros(nvar,length(nstepsLOC))  ;
        SNAP.(fff{iii})(:,1) = VAR_n.(fff{iii})(:,end) ;
    end
end

%
% fff =fieldnames(GAUSSV) ;
% SNAP_GAUSSV = [] ;
%
% for iii=1:length(fff)
%      nvar = size(GAUSSV.(fff{iii}),1) ;
%      SNAP_GAUSSV.(fff{iii}) = zeros(nvar,length(nstepsLOC))  ;
% end


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

NODESV_PROP.ACEL.TYPE = 'Vector' ;  % For printing in GID
NODESV_PROP.ACEL.LEGEND = 'ACEL.' ;  % GID's printing
if DATA.MESH.ndim == 2
    NODESV_PROP.ACEL.COMP = {'a-x','a-y'} ;
else
    NODESV_PROP.ACEL.COMP = {'a-x','a-y','a-z'} ;
end



NODESV_PROP.VELOC.TYPE = 'Vector' ;  % For printing in GID
NODESV_PROP.VELOC.LEGEND = 'VELOC.' ;  % GID's printing
if DATA.MESH.ndim == 2
    NODESV_PROP.VELOC.COMP = {'v-x','v-y'} ;
else
    NODESV_PROP.VELOC.COMP = {'v-x','v-y','v-z'} ;
end




GAUSSV_PROP.PK2STRESS.LEGEND = 'PK2STRESS' ;


if DATA.MESH.ndim  == 2
    GAUSSV_PROP.PK2STRESS.TYPE = 'Matrix' ;
    GAUSSV_PROP.PK2STRESS.COMP ={'S-X','S-Y','S-Z','S-XY'} ;
    switch DATA.typePROBLEM
        case 'pstress'
            GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1,2,4,3] ;
        otherwise
            warning('Stresses in the z-direction have been not computed...')
            GAUSSV_PROP.PK2STRESS.CONVERSION_GID =[1,2,4,3] ;
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
            warning('Stresses in the z-direction have been not computed...')
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


GAUSSV_PROP.VONMISES_CAUCHY_STRESS.LEGEND = 'CAUCHY VONMISES' ;


GAUSSV_PROP.VONMISES_CAUCHY_STRESS.TYPE = 'Scalar' ;


DATA.PRINT.NODESV_PROP = NODESV_PROP  ;
DATA.PRINT.GAUSSV_PROP  =GAUSSV_PROP  ;

DATA_cl.COMMON = DATA ;


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
% GAUSSV.PK2STRESS = zeros(ndofSTRESS,nstoreSTRESS) ;  % 2nd Piola-Kirchhoff stress
% GAUSSV.GLSTRAINS = zeros(ndofSTRESS,nstoreSTRESS) ;  % Green-Lagrange strains
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
