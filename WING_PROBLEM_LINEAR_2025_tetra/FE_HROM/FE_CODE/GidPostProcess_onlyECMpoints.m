function GidPostProcess_onlyECMpoints(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DATA,WHICHMODE);
% Post-processing of results using GID  (vibration modes)
%dbstop('4')
if nargin==0
    load('tmp.mat')
end

% Name of the mesh fileRpo
[DIRinp,NAMEFILEinput ]=fileparts(NameFileMesh) ;

ROOTFOLDER = [DIRinp,filesep,'GIDPOST',filesep] ;
if ~exist(ROOTFOLDER)
    mkdir(ROOTFOLDER) ;
end
NameFile_msh = [ROOTFOLDER,filesep,NAMEFILEinput,'' ,'.msh'] ;

NameFile_res= [ROOTFOLDER,filesep,NAMEFILEinput,'','.res'] ;

% Writing mesh file
%MaterialType = ones(size(CN,1),1) ;
NAMEPROJ = ['MODES_',WHICHMODE,'_'] ; 
DATA = DefaultField(DATA,'NAMEPROJ_BASE',NAMEPROJ) ; 
NAMEPROJ = DATA.NAMEPROJ_BASE ; %  ['MODES_',WHICHMODE,'_'] ; 


DATA = DefaultField(DATA,'MaterialType',ones(size(CN,1),1)) ; 


GidMesh2DFE(NameFile_msh,COOR,CN,NAMEPROJ,DATA.MaterialType,TypeElement);
% Writing results file
DATA = DefaultField(DATA,'nMODESplot',size(MODES,2)) ;
nMODES = DATA.nMODESplot   ;
DATA = DefaultField(DATA,'NODES',1:size(COOR,1)) ;
% if ~isempty(DOFl)
%     MODESplot = zeros(length(DATA.NODES)*size(COOR,2),DATA.nMODESplot) ;
%     
%     MODESplot(DOFl,:) = MODES(:,1:nMODES) ;
% else
    MODESplot = MODES ;
%end
GidResults2DFE_modes(NameFile_res,COOR,CN,TypeElement,MODESplot,posgp,DATA);

% GidResults2DFE(NameFile,COOR,CONNECT,TypeElement,d,strainGLO, stressGLO, ...
%     React,NAME_INPUT_DATA,posgp,DATA)

NAMEFILEOPEN =  [NameFile_res] ;

if isunix 
    switch  (NAMEFILEOPEN(1))
        case '/' 
        otherwise
        NAMEFILEOPEN = [cd,filesep,NAMEFILEOPEN];
    end
end
  
disp('open GID FILE:')
disp(NAMEFILEOPEN)