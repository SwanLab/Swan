function GidPostProcessModesFINT(COOR,CN,TypeElement,MODES,posgp,DATAloc);
% Post-processing of results using GID  (vibration modes)
%dbstop('4')
if nargin==0
    load('tmp1.mat')
end

 


% DATAIN.LABEL_NAME_LOC = 'ErrorSmoothBmat' ; 
% DATAIN.MESSAGEloc = 'ErrorSmoothBmat' ;
%DATA= DefaultField(DATA,'LABEL_NAME_FILE','MOD_fINT') ; 

%    DATAIN.LABEL_NAME_VARIABLE = 'ErrSmoothFint(t.p.)' ; 
%    DATAIN.LABEL_NAME_FILE = 'ErrSmF' ;


% Name of the mesh fileRpo
NameFile_msh = [DATAloc.NameFileMesh,'.msh'] ;
% Name of the results file
NameFile_res= [DATAloc.NameFileMesh,'.res'] ;


 DATAloc.LABEL_NAME_FILE = ' '; 

% Writing mesh file
DATAloc = DefaultField(DATAloc,'MaterialType',ones(size(CN,1),1)) ; 
%MaterialType = ones(size(CN,1),1) ;
GidMesh2DFE(NameFile_msh,COOR,CN,DATAloc.LABEL_NAME_FILE,DATAloc.MaterialType,TypeElement);
% Writing results file
DATAloc = DefaultField(DATAloc,'nMODES',size(MODES,2)) ;
nMODES = DATAloc.nMODES   ;
DATAloc = DefaultField(DATAloc,'NODES',1:size(COOR,1)) ;
% if ~isempty(DOFl)
%     MODES = zeros(length(DATA.NODES)*size(COOR,2),DATA.nMODES) ;
%     
%     MODES(DOFl,:) = MODES(:,1:nMODES) ;
% else
    %MODES = MODES ;
% end
GidResults2DFE_modesFINT(NameFile_res,COOR,CN,TypeElement,MODES,posgp,DATAloc);

cddd = cd ;
%NAMEFILEOPEN =  [cddd,filesep,NameFile_res] ;

NAMEFILEOPEN  = NameFile_res ; 

disp('open GID FILE:')
disp(NAMEFILEOPEN)
 
