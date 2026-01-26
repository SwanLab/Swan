function [NODESV_n,GAUSSV_n,NODESV_PROP,GAUSSV_PROP] = InitHistoricVariables(BND_Dyn,OPERfe,PROPMAT,DATA)


DATA = DefaultField(DATA,'PRINT_GID',[]) ;
DATA = DefaultField(DATA,'STORE',[]) ;
 

NODESV_n.U = BND_Dyn.d0 ;  % Initial displacements
NODESV_n.VEL = BND_Dyn.v0  ; % Nodal velocities (only DOFf)  at time n
NODESV_n.ACEL= zeros(size(NODESV_n.VEL)); % Nodal velocities (only DOFf)  at time n
NODESV_n.Reactions= zeros(size(NODESV_n.VEL)); % Nodal velocities (only DOFf)  at time n

% Properties variables at nodes
% Displacements
% ----------------------------------
NODESV_PROP.U.TYPE = 'Vector' ;
NODESV_PROP.U.PRINT = 1 ;  % GID's printing
NODESV_PROP.U.LEGEND = 'DISPLAC.' ;  % GID's printing
if OPERfe.ndimSP == 2
    NODESV_PROP.U.COMP = {'u-x','u-y'} ;
else
    NODESV_PROP.U.COMP = {'u-x','u-y','u-z'} ;
end

% Reactions
NODESV_PROP.Reactions.TYPE = 'Vector' ;
NODESV_PROP.Reactions.PRINT = 1 ;  % GID's printing
NODESV_PROP.Reactions.STORE = 1 ;  % Store in memory
DATA.PRINT_GID = DefaultField(DATA.PRINT_GID,'REACTIONS',1) ;
 
 
if DATA.PRINT_GID.REACTIONS == 1
    NODESV_PROP.Reactions.PRINT = 1 ;  % GID's printing
else
    NODESV_PROP.Reactions.PRINT =0 ;
     NODESV_PROP.Reactions.ISSPARSE =1 ;
end





NODESV_PROP.Reactions.LEGEND = 'REACT.' ;  % GID's printing
if OPERfe.ndimSP == 2
    NODESV_PROP.Reactions.COMP = {'R-x','R-y'} ;
else
    NODESV_PROP.Reactions.COMP = {'R-x','R-y','R-z'} ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GAUSSV_n.stressST = zeros(OPERfe.ngausT,1);  % stressSTes at time n+1
if OPERfe.ndimSP == 2 && OPERfe.nstrain == 4
    GAUSSV_PROP.stressST.TYPE = 'Matrix' ;
    GAUSSV_PROP.stressST.COMP ={'SIGMA-X','SIGMA-Y','SIGMA-Z','TAU-XY'} ;
    GAUSSV_PROP.stressST.CONVERSION_GID =[1,2,3,4] ;
elseif OPERfe.ndimSP == 3 && OPERfe.nstrain == 6
    GAUSSV_PROP.stressST.TYPE = 'Matrix' ;
    GAUSSV_PROP.stressST.COMP ={'SIGMA-X','SIGMA-Y','SIGMA-Z','TAU-XY','TAU-YZ','TAU-XZ'} ;
    GAUSSV_PROP.stressST.CONVERSION_GID =[1 2 3 6 4 5] ;
end
DATA.PRINT_GID = DefaultField(DATA.PRINT_GID,'STRESSES',1) ;
if DATA.PRINT_GID.STRESSES == 1
    GAUSSV_PROP.stressST.PRINT = 1 ;
else
    GAUSSV_PROP.stressST.PRINT = 0 ;
    
end
DATA.STORE = DefaultField(DATA.STORE,'STRESSES',1) ;
if DATA.STORE.STRESSES == 0
    GAUSSV_PROP.stressST.STORE = 0 ;
end


GAUSSV_PROP.stressST.LEGEND = 'STRESS' ;







switch DATA.TypeImplementation
    case 'J2_plasticity_small_strains'
        GAUSSV_n.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
        GAUSSV_n.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
        GAUSSV_n.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
        
        GAUSSV_PROP.alpha.TYPE = 'Scalar' ;
        GAUSSV_PROP.alpha.PRINT = 1 ;  % GID's printing
        GAUSSV_PROP.alpha.LEGEND = 'INT.VAR.' ;
        
        GAUSSV_n.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
        GAUSSV_n.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
        GAUSSV_n.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
        GAUSSV_n.VonMises = zeros(OPERfe.ngaus,1);; % stressST-like internal variable
        
        
        
        
        GAUSSV_PROP.alpha.TYPE = 'Scalar' ;
        
        DATA.PRINT_GID = DefaultField(DATA.PRINT_GID,'INTERNAL_VARIABLES',1) ;
        if DATA.PRINT_GID.INTERNAL_VARIABLES == 1
            GAUSSV_PROP.alpha.PRINT = 1 ;  % GID's printing
        else
            GAUSSV_PROP.alpha.PRINT = 0 ;
        end
        GAUSSV_PROP.alpha.LEGEND = 'INT.VAR.' ;
        
        
        DATA.STORE = DefaultField(DATA.STORE,'INTERNAL_VARIABLES',1) ;
        if DATA.STORE.INTERNAL_VARIABLES == 0
            GAUSSV_PROP.alpha.STORE = 0 ;
        end
        
        
        GAUSSV_PROP.VonMises.TYPE = 'Scalar' ;
        GAUSSV_PROP.VonMises.PRINT = 1 ;  % GID's printing
        GAUSSV_PROP.VonMises.LEGEND = 'Von Mises' ;
    otherwise
        error('Option not implemented')
end



% -----------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSSV_np1.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
% GAUSSV_np1.stressST = zeros(OPERfe.ngausT,1);  % stressSTes at time n+1
% GAUSSV_np1.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
% GAUSSV_np1.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
% % -----------------------------------------------------------------------
