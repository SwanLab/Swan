function GidPostProcess2D(COOR,CN,TypeElement,d,strainGLO, stressGLO,  React,NAME_INPUT_DATA,posgp,NameFileMesh); 
% Post-processing of results using GID (2D)
if nargin==0
    load('tmp.mat')
end

% Name of the mesh file 
NameFile_msh = ['GIDPOST/',NAME_INPUT_DATA,'_',NameFileMesh(1:end-4),'.msh'] ; 
% Name of the results file 
NameFile_res= ['GIDPOST/',NAME_INPUT_DATA,'_',NameFileMesh(1:end-4),'.res'] ; 

% Writing mesh file
MaterialType = ones(size(CN,1),1) ; 
GidMesh2DFE(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,'TypeElement',TypeElement);
% Writing results file
GidResults2DFE(NameFile_res,COOR,CN,TypeElement,d,strainGLO, stressGLO,  React,NAME_INPUT_DATA,posgp);

cddd = cd ; 
NAMEFILEOPEN =  [cddd,'/',NameFile_res] ; 
disp('open GID FILE:')
disp(NAMEFILEOPEN)