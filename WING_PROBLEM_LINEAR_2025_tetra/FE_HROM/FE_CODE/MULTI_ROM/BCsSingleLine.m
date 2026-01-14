function  [DOFr,dR,G,DOFm,NODES_LINES_ROTATIONS] = BCsSingleLine(NODES_LINES,ilines,DISP,ndim,V,DATAIN,MESH2D,DATA_REFMESH)

if nargin == 0
    load('tmp.mat')
   
end
NODES_LINES_ROTATIONS =  [] ; 

DATAIN = DefaultField(DATAIN,'GLOBAL_DIRICHLET_BOUNDARY_CONDITIONS',0) ; % 28-Jan-2020. BC in global axis


DISPLOC  =DISP.LINE{ilines} ;  % Displacement associated to the interfaces
ISSCC = cellfun(@isempty,DISPLOC) ; 

 
 





ISTR = 2;
if DATAIN.ndimSP == 3
    ISTR=  3 ;
end
ISSCCall = ISSCC ;
ISSCCROT = ISSCC(ISTR+1:end) ; % Rotations DOFs
ISSCC = ISSCC(1:ISTR) ;  % Translational modes

DATAIN = DefaultField(DATAIN,'BoundaryConditionsUnsconstrainedDOFS',0) ;
DOFm = [] ; 
G = [] ; 

if   ~all(ISSCCall==1) &&  ~all(ISSCCall==0) &&  DATAIN.BoundaryConditionsUnsconstrainedDOFS == 1
  
    % New method (19 July 2019)
    
    if DATAIN.GLOBAL_DIRICHLET_BOUNDARY_CONDITIONS  == 1
        error('Option not implemented yet for unconstrainedDOFS')
    end
    
    [DOFr,dR,G,DOFm] = ...
        BCsCoarseScaleGeneral(ISSCCROT,ISSCC,NODES_LINES,ilines,DISPLOC,ndim,V,DATAIN,MESH2D,DATA_REFMESH) ; 
    
    
else
    
    if   DATAIN.GLOBAL_DIRICHLET_BOUNDARY_CONDITIONS  == 1
        
        if all(ISSCC)==0 && any(ISSCC) 
            error('Option not implemented for unconstrained DOFs')
        elseif all(ISSCC) ==1
            DOFr  = []; dR = [] ;  % Nothing is done 
        else
            % Conversion to local coordinates  is necessary (29-Jan-2020)
            % -----------------------------------------------------------
            [DOFr,dR,NODES_LINES_ROTATIONS] =  DirichletBCglobalRVE(NODES_LINES{ilines},ndim,DISPLOC,MESH2D,DATA_REFMESH,V,DATAIN)  ; 
        end
    
        
    else
        % Classic method (first one implemented) 
   [DOFr,dR] =  DirichletBCstandardRVE(NODES_LINES,ilines,ndim,DISPLOC)  ; 

    end
    
    
    
end



