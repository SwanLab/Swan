function [VintfDEF,COORbnd,CNbREN,b,MintfCHOL] = DefPartInterfaceModes(MESH,PhiRB,Vintf,Mintf,DATA,nDEF_basic,DATAoffline) 
% Deformational component of fictititious interface modes
% 7-Apr-2024, Sunday, Molinos Marfagones, Cartagena 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx

% Vintf
% 1. Coordinate boundary nodes
disp(['Deformational component of the interface modes'])
BoundaryNodes = MESH.faceNODESall ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;
b = MESH.faceDOFSall;
VintfRB = PhiRB(b,:) ;
VintfDEF = SprojDEF_operator(VintfRB,Mintf,Vintf) ;
DATAloch.TOL = 1e-3;


%DATAoffline = DefaultField(DATAoffline,'USE_CHOLESKY_DECOMPOSITION',1) ; % JAHO 22-Apr-2024

if DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
[ VintfDEF,Sss,Vss,MintfCHOL] = WSVDT( VintfDEF,Mintf,DATAloch) ;
else
    [ VintfDEF,Sss,Vss,MintfCHOL] = WSVDT( VintfDEF,[],DATAloch) ;

end
disp(['Number of deformational interface modes = ', num2str(size(VintfDEF,2)),'  (see GID plot below)'])




CNb = cell(length(MESH.BNDinterface),1) ;
for iface  = 1:length(MESH.BNDinterface)
    CNb{iface} = MESH.BNDinterface(iface).CNb ;
end
CNb = cell2mat(CNb(:)) ;
CNbREN  =  RenumberConnectivities( CNb,1:length(BoundaryNodes) );
% ------------------------------------------------------------------------------------------------------------------------
NameLoc =     'DispIntf_DEF' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [];
disp(['Deformed shapes of deformational interface modes (GID) --> '])
GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[VintfDEF],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC) ;
disp('-------------------------------------')