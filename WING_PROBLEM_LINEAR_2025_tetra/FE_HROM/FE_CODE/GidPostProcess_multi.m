function GidPostProcess_multi(COOR,CN,TypeElement,d,strainGLO, stressGLO, ...
    React,NAME_INPUT_DATA,posgp,NameFileMesh,MaterialType,DATA,NODES,RESIDUAL,DATAINM);
% Post-processing of results using GID
%dbstop('5')
if nargin==0
    load('tmp2.mat')
end

% Name of the mesh file
%if ~isempty(NameFileMesh)
%[dummy1 NameFileMeshHERE dummy2]= fileparts(NameFileMesh) ;
%else
NameFileMeshHERE = '' ;
%end

DATA = DefaultField(DATA,'POST_LabelFileGid','') ;

NameFile_msh = ['GIDPOST',filesep,NameFileMeshHERE,NAME_INPUT_DATA,'_',DATA.POST_LabelFileGid ,'.msh'] ;
% Name of the results file
NameFile_res= ['GIDPOST',filesep,NameFileMeshHERE,NAME_INPUT_DATA,'_',DATA.POST_LabelFileGid,'.res'] ;


IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,TypeElement,DATA.NAMEMESH);
% Writing results file
% Nodal vectors
% -------------------
NODAL_VECTOR(1).NAME = 'DISPLACEMENTS' ;
NODAL_VECTOR(1).COMP = {'X-DISP','Y-DISP','Z-DISP'};
NODAL_VECTOR(1).VAR = d ;
NODAL_VECTOR(1).NODES = 1:size(COOR,1) ;
%
%
if ~isempty(RESIDUAL.ALL)
    if DATAINM.FACTOR_MULTIPLY_RESIDUAL_FORCES > 1
        APPLEG = ['(x',num2str(DATAINM.FACTOR_MULTIPLY_RESIDUAL_FORCES),')'] ; 
    end
    
    
    NODAL_VECTOR(2).NAME = ['RES.FORCES',APPLEG] ;
    NODAL_VECTOR(2).COMP = {'X-RES','Y-RES','Z-RES'};
    NODAL_VECTOR(2).VAR = RESIDUAL.ALL ;
    NODAL_VECTOR(2).NODES = NODES{2}' ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scalar vectors: forces and moments
DATA = DefaultField(DATA,'Forces1D',[]) ;
NODAL_SCALAR = [] ; iacum = 1 ;
if ~isempty(DATA.Forces1D)
     if DATAINM.FACTOR_MULTIPLY_MOMENTS_DIAGRAM > 1
        APPLEG = ['(x',num2str(DATAINM.FACTOR_MULTIPLY_MOMENTS_DIAGRAM),')'] ; 
    end
    
    fff = fieldnames(DATA.Forces1D) ;
    for ifield = 1:length(fff)
        FORCE = DATA.Forces1D.(fff{ifield}) ;
        if ~isempty(FORCE)
            NODAL_SCALAR(iacum).NAME = [fff{ifield},APPLEG ];
            NODAL_SCALAR(iacum).VAR = FORCE ;
            NODAL_SCALAR(iacum).NODES = NODES{1} ;
            iacum = iacum + 1;
        end
    end
end
% scalar defined on GAUSS POINTS, mesh 1
imesh = 1 ;
MESH(imesh).TypeElement = TypeElement{imesh} ;
MESH(imesh).posgp = posgp{imesh} ;
MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
MESH(imesh).GAUSS_SCALAR(1).NAME = 'Max.Von_Mises' ;
MESH(imesh).GAUSS_SCALAR(1).VAR = DATA.MAXstressVONMISES ;
MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;

if ~isempty(RESIDUAL.MAX_NORM)
     if DATAINM.FACTOR_MULTIPLY_RESIDUAL_FORCES > 1
        APPLEG = ['(x',num2str(DATAINM.FACTOR_MULTIPLY_RESIDUAL_FORCES),')'] ; 
    end
MESH(imesh).GAUSS_SCALAR(2).NAME = ['RESIDF_Max_Norm',APPLEG] ;
MESH(imesh).GAUSS_SCALAR(2).VAR = RESIDUAL.MAX_NORM ;
MESH(imesh).GAUSS_SCALAR(2).ELEMENTS = IND_ELEM_MESHES{imesh} ;
MESH(imesh).GAUSS_SCALAR(3).NAME = ['RESIDF_Norm_Avg',APPLEG] ;
MESH(imesh).GAUSS_SCALAR(3).VAR = RESIDUAL.NORM_AVERAGE ;
MESH(imesh).GAUSS_SCALAR(3).ELEMENTS = IND_ELEM_MESHES{imesh} ;
end



imesh = 2;
MESH(imesh).TypeElement = TypeElement{imesh} ;
MESH(imesh).posgp = posgp{imesh} ;
MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
MESH(imesh).GAUSS_SCALAR(1).NAME = 'Von_Mises' ;
MESH(imesh).GAUSS_SCALAR(1).VAR = DATA.stressVONMISES ;
MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
%
MESH(imesh).GAUSS_MATRIX(1).NAME = 'STRESS' ;
MESH(imesh).GAUSS_MATRIX(1).VAR = stressGLO ;
MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
MESH(imesh).GAUSS_MATRIX(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;

nelemE(1) = size(CN{1},2) ;
nelemE(2) = size(CN{2},2) ;
ndim = size(COOR,2) ;

GidResults2DFE_multi(NameFile_res,ndim,TypeElement,NODAL_VECTOR,NODAL_SCALAR,MESH,posgp);

cddd = cd ;
NAMEFILEOPEN =  [cddd,filesep,NameFile_res] ;
NAMEMESH = [cddd,filesep,NameFile_msh] ;

%  if DATA.COMPRESS_GID == 1
%     disp('Compressing GID FILE ...')
%     TEXT = ['gid -PostResultsToBinary  ',NAMEFILEOPEN, ' ',NAMEFILEOPEN] ;
%     unix(TEXT) ;
%     delete(NameFile_msh)
%     disp('Done')
% end


disp('open GID FILE:')
clipboard('copy',NAMEFILEOPEN)
disp(NAMEFILEOPEN)

DATA = DefaultField(DATA,'OPEN_GID',0) ;
if DATA.OPEN_GID ==1
    TTTT = ['gidpost ',NAMEFILEOPEN] ;
    unix(TTTT);
end