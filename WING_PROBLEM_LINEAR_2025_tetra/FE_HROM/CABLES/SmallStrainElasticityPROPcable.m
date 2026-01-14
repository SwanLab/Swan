function [MATPRO] = SmallStrainElasticityPROPcable(MESH,PROPMAT)

if nargin == 0
    load('tmp.mat')
end

ndim = size(MESH.COOR,2)  ;

MESH = DefaultField(MESH,'nstrain',[]) ;
% if isempty(MESH.nstrain)
%     if ndim==2
%         nstrain = 3;
%     else
%         nstrain = 6 ;
%         typePROBLEM ='3D' ;
%     end
% else
nstrain = MESH.nstrain ;
% end

nelem = size(MESH.MaterialType,1) ;
MATPRO.celasglo = zeros(nelem,1) ;  % Global array of elasticity matrices
MATPRO.dens = zeros(nelem,1) ;
%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
    % celas3D =PROPMAT(imat).ElasticityMatrix ; %
    % INVcelas3D = inv(celas3D) ;
    ELEMS = find(MESH.MaterialType == imat) ;
    MATPRO.celasglo(ELEMS) = PROPMAT(imat).EA ;
    MATPRO.dens(ELEMS) = PROPMAT(imat).dens ;
end
