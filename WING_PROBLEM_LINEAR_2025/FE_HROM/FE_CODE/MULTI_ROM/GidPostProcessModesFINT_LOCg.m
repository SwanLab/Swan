function DATA = GidPostProcessModesFINT_LOCg(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DATA);
% Post-processing of results using GID  (vibration modes)
%dbstop('4')
if nargin==0
    load('tmp1.mat')
end

[aaa,bbb] = fileparts(NameFileMesh)   ; 

 
NAMEF = [aaa,filesep,'GIDPOST'] ; 

if ~exist(NAMEF)
    mkdir(NAMEF)
end


% DATAIN.LABEL_NAME_LOC = 'ErrorSmoothBmat' ; 
% DATAIN.MESSAGEloc = 'ErrorSmoothBmat' ;
DATA= DefaultField(DATA,'LABEL_NAME_FILE','MOD_fINT') ; 

%    DATAIN.LABEL_NAME_VARIABLE = 'ErrSmoothFint(t.p.)' ; 
%    DATAIN.LABEL_NAME_FILE = 'ErrSmF' ;


% Name of the mesh fileRpo
NameFile_msh = [NAMEF,filesep,DATA.LABEL_NAME_FILE,'_',bbb(1:end-4),'.msh'] ;
% Name of the results file
NameFile_res= [NAMEF,filesep,DATA.LABEL_NAME_FILE,'_',bbb(1:end-4),'.res'] ;


 

% Writing mesh file
DATA = DefaultField(DATA,'MaterialType',ones(size(CN,1),1)) ; 
%MaterialType = ones(size(CN,1),1) ;
GidMesh2DFE(NameFile_msh,COOR,CN,DATA.LABEL_NAME_FILE,DATA.MaterialType,TypeElement);
% Writing results file
DATA = DefaultField(DATA,'nMODES',size(MODES,2)) ;
nMODES = DATA.nMODES   ;
DATA = DefaultField(DATA,'NODES',1:size(COOR,1)) ;
% if ~isempty(DOFl)
%     MODES = zeros(length(DATA.NODES)*size(COOR,2),DATA.nMODES) ;
%     
%     MODES(DOFl,:) = MODES(:,1:nMODES) ;
% else
    %MODES = MODES ;
% end
GidResults2DFE_modesFINT(NameFile_res,COOR,CN,TypeElement,MODES,posgp,DATA);

cddd = cd ;
%NAMEFILEOPEN =  [cddd,filesep,NameFile_res] ;

NAMEFILEOPEN  = NameFile_res ; 

disp('open GID FILE:')
disp(NAMEFILEOPEN)

DATA = DefaultField(DATA,'MSGPRINT',{}) ;
DATA.MSGPRINT{end+1} = '-----------------------------------------' ; 
DATA.MSGPRINT{end+1} = [DATA.LABEL_NAME_FILE,' --> OPEN GID FILE'] ; 
DATA.MSGPRINT{end+1} =  NAMEFILEOPEN; 
DATA.MSGPRINT{end+1} = '-----------------------------------------' ; 
