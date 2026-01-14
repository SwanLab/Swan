function MESH = NormalsAndTangents(MESH)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/FEoperators_1dom.m

if nargin == 0
    load('tmp.mat')
end
 
  

% Connectivities faces 
 
%Indexes_faces_bnd_element = cell(1,length(MESH.NODES_FACES)) ;
NORMALv = cell(1,length(MESH.CNb)) ; 
TANGENTv = cell(1,length(MESH.CNb)) ; 
for iface = 1:length(MESH.CNb)     
    [NORMALv{iface},TANGENTv{iface}] = NormalsBoundaryLocal(MESH.COOR,MESH.CNb{iface} )  ;     
end

MESH.NormalBoundaryElementsFace = NORMALv ; 
MESH.TangentBoundaryElementsFace = TANGENTv ; 

% % Moved inside on 22-March-2021
% %***********************************
% 
% DATA = DefaultField(DATA,'RenumberElementsForEficiency',1) ;  
% if  DATA.RenumberElementsForEficiency == 1
%     disp('Renumering elements')
%     [~,IndicesRenumberingElements]  = sort(MESH.CN(:,1)) ;
%     MESH.CN = MESH.CN(IndicesRenumberingElements,:) ;
%   %  MATPRO.celasglo = MATPRO.celasglo(:,:,IndicesRenumberingElements) ;
%   %  MATPRO.dens = MATPRO.dens(IndicesRenumberingElements) ;
%     MESH.MaterialType = MESH.MaterialType(IndicesRenumberingElements) ;
% else
%     IndicesRenumberingElements = [] ;
% end
% MESH.IndicesRenumberingElements =IndicesRenumberingElements ; 
% 
