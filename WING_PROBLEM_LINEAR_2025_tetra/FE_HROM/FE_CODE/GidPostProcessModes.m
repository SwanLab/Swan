
function GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DATA,DOFl);
% Post-processing of results using GID  (vibration modes)
%dbstop('4')
if nargin==0
    load('tmp.mat')
end

[aaa,bbb] = fileparts(NameFileMesh)   ; 

if ~isempty(aaa)
NAMEF = [aaa,filesep,'GIDPOST'] ; 
else
    NAMEF = ['GIDPOST'] ; 

end

if ~exist(NAMEF)
    mkdir(NAMEF)
end

% Name of the mesh fileRpo
NameFile_msh = [NAMEF,filesep,'MODES','_',bbb(1:end-4),'.msh'] ;
% Name of the results file
NameFile_res= [NAMEF,filesep,'MODES','_',bbb(1:end-4),'.res'] ;

% Writing mesh file
DATA = DefaultField(DATA,'MaterialType',ones(size(CN,1),1)) ; 
%MaterialType = ones(size(CN,1),1) ;
GidMesh2DFE(NameFile_msh,COOR,CN,'MODES',DATA.MaterialType,TypeElement);
% Writing results file
DATA = DefaultField(DATA,'nMODESplot',size(MODES,2)) ;
nMODES = DATA.nMODESplot   ;
DATA = DefaultField(DATA,'NODES',1:size(COOR,1)) ;
if ~isempty(DOFl)
    MODESplot = zeros(length(DATA.NODES)*size(COOR,2),DATA.nMODESplot) ;
    
    MODESplot(DOFl,:) = MODES(:,1:nMODES) ;
else
    MODESplot = MODES ;
end
GidResults2DFE_modes(NameFile_res,COOR,CN,TypeElement,MODESplot,posgp,DATA);

cddd = cd ;
%NAMEFILEOPEN =  [cddd,filesep,NameFile_res] ;

NAMEFILEOPEN =  [NameFile_res] ;

disp('open GID FILE:')
disp(NAMEFILEOPEN)