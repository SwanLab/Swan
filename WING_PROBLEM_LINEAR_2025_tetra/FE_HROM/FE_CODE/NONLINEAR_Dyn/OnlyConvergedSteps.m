function [GAUSS_SNAP,NODES_SNAP,DATA] = OnlyConvergedSteps(GAUSS_SNAP,NODES_SNAP,istep_conv,DATA)

if nargin == 0
    load('tmp1.mat')
end


NAMEFIELDS = fieldnames(NODES_SNAP) ; 
for inames = 1:length(NAMEFIELDS)
    nameloc = NAMEFIELDS{inames} ; 
    NODES_SNAP.(nameloc) = NODES_SNAP.(nameloc)(:,1:istep_conv) ; 
end
 

NAMEFIELDS = fieldnames(GAUSS_SNAP) ; 
for inames = 1:length(NAMEFIELDS)
    nameloc = NAMEFIELDS{inames} ; 
    GAUSS_SNAP.(nameloc) = GAUSS_SNAP.(nameloc)(:,1:istep_conv)  ; 
end

%DATA.TIME_DISCRETIZATION = DATA.TIME_DISCRETIZATION(1:istep_conv) ; 
 