function [DATAIN] = PrintingGid_ECMpoints_LARGE(MESH,DATA_GENGAUSS,HISTORY_ITERATIONS) ;

if nargin == 0
    load('tmp1.mat')
end


ELEMS_ini = HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{1} ; % initial set of elements
ELEMS_FIN = HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{end} ; % Final set of elements
wINI = HISTORY_ITERATIONS.WEIGHTS{1} ;
wFIN = HISTORY_ITERATIONS.WEIGHTS{end} ;

[elemsREMAAIN,indOLDremain,indNEWremain ]= intersect(ELEMS_ini,ELEMS_FIN) ;
[elemsNEW,indNEWnew]= setdiff(ELEMS_FIN,ELEMS_ini) ;
indOLDdissapear = 1:length(ELEMS_ini) ;
indOLDdissapear(indOLDremain)= [] ;

HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS = {} ;

HISTORY_ITERATIONS.WEIGHTS = {} ;
imesh = 1 ;
NAMEMESH{imesh} = 'INI_removed' ;
HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{imesh} = ELEMS_ini(indOLDdissapear);
HISTORY_ITERATIONS.WEIGHTS{imesh} =wINI(indOLDdissapear) ;
imesh = 2 ;
NAMEMESH{imesh} = 'INI_remain' ;
HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{imesh} = elemsREMAAIN;
HISTORY_ITERATIONS.WEIGHTS{imesh} =wINI(indOLDremain) ;
imesh = 3 ;
NAMEMESH{imesh} = 'FINAL_remain' ;
HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{imesh} = elemsREMAAIN;
HISTORY_ITERATIONS.WEIGHTS{imesh} =wFIN(indNEWremain) ;
imesh = 4 ;
NAMEMESH{imesh} = 'FINAL_new' ;
HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{imesh} = elemsNEW;
HISTORY_ITERATIONS.WEIGHTS{imesh} =wFIN(indNEWnew) ;

imesh = 5 ;
NAMEMESH{imesh} = 'FULL' ;




DATA= [] ;


% 
% % Name of the mesh fileRpo
% [DIRinp,NAMEFILEinput ]=fileparts(DATAIN.NAME_WS_MODES) ;
% DATAIN = DefaultField(DATAIN,'LABEL_NAME_PROJECT','') ;
% NAMEFILEinput=  ['ECM_ELEMENTS',DATAIN.LABEL_NAME_PROJECT] ;
% 
% ROOTFOLDER = [DIRinp,filesep,'GIDPOST',filesep] ;
% disp(['Writing ECM mesh  in folder: ',ROOTFOLDER]);
% if ~exist(ROOTFOLDER)
%     mkdir(ROOTFOLDER) ;
% end
NameFile_msh = [DATA_GENGAUSS.NameFileMesh_CECM,'' ,'.msh'] ;
NameFile_res = [DATA_GENGAUSS.NameFileMesh_CECM,'' ,'.res'] ;


% Writing mesh file
%MaterialType = ones(size(CN,1),1) ;
CNred = cell(size(HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS)) ;




MatRED =CNred ;
for i = 1:length(CNred)
    CNred{i}  = MESH.CN(HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{i},:) ;
    MatRED{i} = i*ones(size(MESH.MaterialType(HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{i})))  ;
end




COORprint = MESH.COOR ;
%CNprint = {CNred,CNref} ;
CNprint = CNred ; CNprint{end+1} = MESH.CN;
NAME_INPUT_DATA =' '  ;
MaterialTypegloPRINT = MatRED ; MaterialTypegloPRINT{end+1} = MESH.MaterialType;
% MaterialTypegloPRINT = {MatRED,DATA_REFMESH.MaterialType} ;
TypeElementPRINT = cell(size(CNprint)) ; % {DATA_REFMESH.TypeElement,DATA_REFMESH.TypeElement} ;


for  i = 1:length(NAMEMESH)
    % NAMEMESH{i} = ['ITER',num2str(i)] ;
    TypeElementPRINT{i} = MESH.TypeElement  ;
end
% NAMEMESH{end} = ['FULL'] ;
% TypeElementPRINT{end} = DATA_REFMESH.TypeElement ;
% %end





IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COORprint,CNprint,NAME_INPUT_DATA,MaterialTypegloPRINT,...
    TypeElementPRINT,NAMEMESH);
%    WEIGHTS = HYPERREDUCED_VARIABLES.WdomRED ;




DATA.NAMEMESH = NAMEMESH ;
%     DATA.WEIGHTS = WEIGHTS;
posgp = cell(size(TypeElementPRINT));
DATAINloc = [] ;
WEIGHTS = cell(size(posgp))  ;
for i = 1:length(WEIGHTS)-1
    WEIGHTS{i} = HISTORY_ITERATIONS.WEIGHTS{i}/sum(HISTORY_ITERATIONS.WEIGHTS{i})*100 ;
end

[NameFileRES ]= GidPostProcess_LARGE_weights(COORprint,CNprint,TypeElementPRINT, ...
    posgp,NameFile_res,MaterialTypegloPRINT,IND_ELEM_MESHES,...
    MESH,DATA,WEIGHTS);




disp('-----------------------------')
disp(['To see reduced set of integration points, open']);
disp(NameFile_msh);
disp(['--------------------- '])
% 
% DATAIN = DefaultField(DATAIN,'MSGPRINT',{});
% DATAIN.MSGPRINT{end+1} = '---------------------------------------------------' ;
% DATAIN.MSGPRINT{end+1} = ['To see reduced set of integration points, open'] ;
% DATAIN.MSGPRINT{end+1} =NameFile_msh ;


