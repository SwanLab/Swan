function Vall_deform = DefComponentInterfaceModes(MESH,PhiRB,Vall,Mintf,INFOPLOTMESHBND,DATA)
% DEFORMATIONAL PART OF THE DISPLACEMENT INTERFACE MODES Vall
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx
% JAHO, 2-April-2023
% ----------------------------------
if nargin == 0
    load('tmp.mat')
end
%  DEFORMATIONAL PART OF THE DISPLACEMENT INTERFACE MODES Vall
PhiRB_f = PhiRB(MESH.faceDOFSall,:) ;  % RIGID-BODY MODES AT THE INTERFACE
% Let us begin by orthogonalizing Vall
Vall_orth = WSVDT(Vall,Mintf) ;
DATALOC.TOLrelSVD = 1e-4; 
[Vall_deform,SS ]=   PprojDEF_operator(PhiRB_f,Mintf,Vall_orth,DATALOC);
% All the singular values associated to deformational modes should equal to
% 1. The remaining ones should be significantly less than one
nmodesORIG = size(Vall,2)  ;  % Number of interface modes
nmodesDEF = size(Vall_deform,2) ;  % Number of deformational modes identified by projection
TOLsingv = 0.9 ;
RBMindex = find(SS<TOLsingv) ;  % INDEXES non-one singular values
ndim = size(MESH.COOR,2) ;  % Number of spatial dimensions
if ndim == 3
    if (nmodesORIG-nmodesDEF)+ length(RBMindex) ~=6
        error('Rigid-body modes have not been properly identified')
    end
else
    if (nmodesORIG-nmodesDEF)+ length(RBMindex) ~=3
        error('Rigid-body modes have not been properly identified')
    end
end
DEFINDEX =setdiff(1:nmodesDEF,RBMindex) ;
Vall_deform = Vall_deform(:,DEFINDEX)  ;


if ~isempty(INFOPLOTMESHBND)
    NameLoc =     'DefPartIntfModes' ;
    NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
    NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
    NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
    NameFile_res = [NAME_MODES_DISP,'.res'] ;
    DATALOC = [];
    
    Vall_deformPLOT = Vall_deform ; 
    nv = sqrt(sum(Vall_deform.^2,1)) ;
    Vall_deformPLOT = bsxfun(@times,Vall_deformPLOT',1./nv')' ; 
    
    GidPostProcessModesDOML(INFOPLOTMESHBND.COOR, INFOPLOTMESHBND.CN ,INFOPLOTMESHBND.TypeElement,Vall_deformPLOT,[],NameFileMesh,NameFile_res,[],DATALOC) ;
%     INFOPLOTMESHBND.COOR = COORbnd;
%     INFOPLOTMESHBND.CN = CNbREN;
%     INFOPLOTMESHBND.TypeElement = MESH.TypeElementB;
%     INFOPLOTMESHBND.posgp = [];
%     INFOPLOTMESHBND.NameFileMesh = NameFileMesh;
%     INFOPLOTMESHBND.NameFile_res = NameFileMesh;
    
    
end