function [NODES_SNAP,GAUSS_SNAP] = InitSnapShotMatrices(OPERfe,DynDATA,DATA,NODESV_n,GAUSSV_n,NODESV_PROP,GAUSSV_PROP) 

% Initialization SNAPSHOT matrices (Gauss variabes, nodal variables)
if nargin == 0
    load('tmp1.mat')
end


nstepsLOC = length(DATA.STEPS_TO_STORE) ; 


% NODAL VARIABLES to be stored in memory 
% ---------------------------------------
NODES_SNAP = [] ; 
VARIABLES = fieldnames(NODESV_PROP) ; 
for ivar = 1:length(VARIABLES)
    NAMEVAR = VARIABLES{ivar} ; 
    VAR_0 = NODESV_n.(NAMEVAR) ; 
    ncomp = length(VAR_0) ; 
    CURRENT_FIELD =  NODESV_PROP.(NAMEVAR) ; 
    CURRENT_FIELD = DefaultField(CURRENT_FIELD,'ISSPARSE',0) ; 
    if CURRENT_FIELD.ISSPARSE == 0
        NODES_SNAP.(NAMEVAR) = zeros(ncomp,nstepsLOC) ;
        NODES_SNAP.(NAMEVAR)(:,1) = VAR_0 ;
    else
        NODES_SNAP.(NAMEVAR) = sparse(ncomp,nstepsLOC) ;
        NODES_SNAP.(NAMEVAR)(:,1) = sparse(VAR_0) ;
    end
end

% GAUSS VARIABLES to be stored in memory 
% ---------------------------------------
GAUSS_SNAP = [] ; 
VARIABLES = fieldnames(GAUSSV_PROP) ; 
for ivar = 1:length(VARIABLES)
    NAMEVAR = VARIABLES{ivar} ; 
    VAR_0 = GAUSSV_n.(NAMEVAR) ; 
    ncomp = length(VAR_0) ; 
    CURRENT_FIELD =  GAUSSV_PROP.(NAMEVAR) ; 
    CURRENT_FIELD = DefaultField(CURRENT_FIELD,'STORE',1) ;
    if CURRENT_FIELD.STORE == 1
    GAUSS_SNAP.(NAMEVAR) = zeros(ncomp,nstepsLOC) ; 
    GAUSS_SNAP.(NAMEVAR)(:,1) = VAR_0 ; 
    else
        GAUSS_SNAP.(NAMEVAR)  = [] ; 
    end
end

 