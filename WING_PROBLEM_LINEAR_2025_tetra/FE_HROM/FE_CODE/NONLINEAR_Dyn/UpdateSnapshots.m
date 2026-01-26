function [GAUSS_SNAP,NODES_SNAP] = UpdateSnapshots(GAUSSV_np1,NODESV_np1,istep,...
    GAUSS_SNAP,NODES_SNAP)

if nargin == 0
    load('tmp1.mat')
end


NAMEFIELDS = fieldnames(NODES_SNAP) ; 
for inames = 1:length(NAMEFIELDS)
    nameloc = NAMEFIELDS{inames} ; 
    if ~issparse(NODES_SNAP.(nameloc))
    NODES_SNAP.(nameloc)(:,istep) = NODESV_np1.(nameloc) ; 
    else
        NODES_SNAP.(nameloc)(:,istep) = sparse(NODESV_np1.(nameloc)) ; 
    end
end
 

NAMEFIELDS = fieldnames(GAUSS_SNAP) ; 
for inames = 1:length(NAMEFIELDS)
    nameloc = NAMEFIELDS{inames} ; 
    if ~isempty(GAUSS_SNAP.(nameloc))
    GAUSS_SNAP.(nameloc)(:,istep) = GAUSSV_np1.(nameloc) ; 
    end
end
 