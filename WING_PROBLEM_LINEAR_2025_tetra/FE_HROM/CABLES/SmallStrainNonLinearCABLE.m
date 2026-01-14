function [MATPRO] = SmallStrainNonLinearCABLE(MESH,PROPMAT)

if nargin == 0
    load('tmp.mat')
end

% ndim = size(MESH.COOR,2)  ;
% 
% MESH = DefaultField(MESH,'nstrain',[]) ;
% % if isempty(MESH.nstrain)
% %     if ndim==2
% %         nstrain = 3;
% %     else
% %         nstrain = 6 ;
% %         typePROBLEM ='3D' ;
% %     end
% % else
% nstrain = MESH.nstrain ;
% end

nelem = size(MESH.MaterialType,1) ;
MATPRO.dens = zeros(nelem,1) ;
MATPRO.AREA = zeros(nelem,1) ;
MATPRO.InitialTension = zeros(nelem,1) ;

%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
    % celas3D =PROPMAT(imat).ElasticityMatrix ; %
    % INVcelas3D = inv(celas3D) ;
    ELEMS = find(MESH.MaterialType == imat) ;
    MATPRO.dens(ELEMS) = PROPMAT(imat).density_linear ;
   MATPRO.AREA(ELEMS)  =  PROPMAT(imat).AREA ; % Not a property material, but included here for simplicity ---no variations of area can be handled 
   % unless specifically coded 
   MATPRO.InitialTension(ELEMS) = PROPMAT(imat).InputsCONSTEQ.InitialStress*PROPMAT(imat).AREA ;
    
end
 
    MATPRO.PROPMAT = PROPMAT; 
 
