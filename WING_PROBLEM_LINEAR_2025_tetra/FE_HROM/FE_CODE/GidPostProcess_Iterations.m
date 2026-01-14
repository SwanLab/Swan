function GidPostProcess_Iterations(COOR,CN,TypeElement,DISP,posgp,NameFileMesh,DATA);
% Post-processing of results using GID  (vibration DISP)
%dbstop('4')
if nargin==0
    load('tmp3.mat')
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
NameFile_msh = [NAMEF,filesep,'','_',bbb(1:end-4),'.msh'] ;
% Name of the results file
NameFile_res= [NAMEF,filesep,'','_',bbb(1:end-4),'.res'] ;

% Writing mesh file
DATA = DefaultField(DATA,'MaterialType',ones(size(CN,1),1)) ; 
%MaterialType = ones(size(CN,1),1) ;

DISP_ini = DISP(:,1) ; 
niter = size(DISP,2)-1; 
DISP_NONCONVERGED = DISP(:,2:end) -repmat(DISP_ini,1,niter) ;  


ndim = size(COOR,2)  ; 
DISP_ini = reshape(DISP_ini,ndim,[]) ; 
COOR = COOR + DISP_ini' ; 


GidMesh2DFE(NameFile_msh,COOR,CN,'DISP',DATA.MaterialType,TypeElement);
% Writing results file
% DATA = DefaultField(DATA,'nDISPplot',size(DISP,2)) ;
% nDISP = DATA.nDISPplot   ;
% DATA = DefaultField(DATA,'NODES',1:size(COOR,1)) ;
% if ~isempty(DOFl)
%     DISPplot = zeros(length(DATA.NODES)*size(COOR,2),DATA.nDISPplot) ;
%     
%     DISPplot(DOFl,:) = DISP(:,1:nDISP) ;
% else
%     DISPplot = DISP ;
% end
GidResults2DFE_iter(NameFile_res,COOR,CN,TypeElement,DISP_NONCONVERGED,posgp,DATA);

cddd = cd ;
%NAMEFILEOPEN =  [cddd,filesep,NameFile_res] ;

NAMEFILEOPEN =  [cddd,filesep,NameFile_res] ;

disp('open GID FILE:')
disp(NAMEFILEOPEN)