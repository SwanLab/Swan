function GidPostProcessDynamic(COOR,CN,TypeElement,DISP,NAME_INPUT_DATA,posgp,NameFileMesh,t); 
% Post-processing of results using GID  
%dbstop('4')
if nargin==0
    load('tmp.mat')
end

% Name of the mesh file 
NameFile_msh = ['GIDPOST/',NAME_INPUT_DATA,'_',NameFileMesh(1:end-4),'.msh'] ; 
% Name of the results file 
NameFile_res= ['GIDPOST/',NAME_INPUT_DATA,'_',NameFileMesh(1:end-4),'.res'] ; 

% Writing mesh file
MaterialType = ones(size(CN,1),1) ; 
GidMesh2DFE(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,TypeElement);
% Writing results file
GidResults2DFEDyn(NameFile_res,COOR,CN,TypeElement,DISP,[], [],  [],NAME_INPUT_DATA,posgp,t);

cddd = cd ; 
NAMEFILEOPEN =  [cddd,'/',NameFile_res] ; 
disp('open GID FILE:')
disp(NAMEFILEOPEN)