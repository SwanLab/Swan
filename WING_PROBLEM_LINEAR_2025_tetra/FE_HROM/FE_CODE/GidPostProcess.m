function GidPostProcess(COOR,CN,TypeElement,d,strainGLO, stressGLO, ...
    React,NAME_INPUT_DATA,posgp,NameFileMesh,MaterialType,DATA,Fnodes);
% Post-processing of results using GID. See GidPostProcess.pdf 
%dbstop('5')
if nargin==0
    load('tmp2.mat')     
end


DATA = DefaultField(DATA,'PostProcessWithNOSLICES',1) ;
DATA = DefaultField(DATA,'MATERIAL',[]) ;
DATA = DefaultField(DATA,'MakeMeshByRepetition',[]) ;
DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'nDOM',1) ;
DATA = DefaultField(DATA,'MaterialTypeNoSlices',[]) ;

if isempty(DATA.MakeMeshByRepetition.nDOM)
    DATA.MakeMeshByRepetition.nDOM  =1  ;
end

if DATA.PostProcessWithNOSLICES == 1 
    
    if ~isempty(DATA.MATERIAL) && (DATA.MakeMeshByRepetition.nDOM(1)>1 ) &&  isempty(DATA.MaterialTypeNoSlices)
        nmat = length(DATA.MATERIAL.PLY) ;
        ndom = prod(DATA.MakeMeshByRepetition.nDOM); %(1) ;
        TOTnmat = nmat*ndom ;
        MAT_TypeORIG = (1:TOTnmat) ;
        MAT_TypeNOSLICES = repmat((1:nmat)',ndom,1) ;
        NewMaterialType = zeros(size(MaterialType));
        
        for imat = 1:length(MAT_TypeORIG)
            INDLOC =find(MaterialType==imat) ;
            NewMaterialType(INDLOC) = MAT_TypeNOSLICES(imat) ;
        end
        MaterialType =  NewMaterialType ;
    elseif ~isempty(DATA.MaterialTypeNoSlices)
        % New version --- Meshes with different reference slices , see
        % GidPostProcess.xoj
        MaterialType = DATA.MaterialTypeNoSlices; 
        
    end
end

[NameFile_msh,NameFile_res] = NameFilePostProcessGid(DATA,NAME_INPUT_DATA) ; 




% Writing mesh file
if isempty(MaterialType)
    MaterialType = ones(size(CN,1),1) ;
end
GidMesh2DFE(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,TypeElement);
% Writing results file
DATA.Fnodes = Fnodes ; 
GidResults2DFE(NameFile_res,COOR,CN,TypeElement,d,strainGLO, stressGLO,  React,NAME_INPUT_DATA,posgp,DATA);
 
NAMEFILEOPEN =  [NameFile_res] ;

if isunix 
    switch  (NAMEFILEOPEN(1))
        case '/' 
        otherwise
        NAMEFILEOPEN = [cd,filesep,NAMEFILEOPEN];
    end
end
  

 
%  if DATA.COMPRESS_GID == 1
%     disp('Compressing GID FILE ...')
%     TEXT = ['gid -PostResultsToBinary  ',NAMEFILEOPEN, ' ',NAMEFILEOPEN] ;
%     unix(TEXT) ;
%     delete(NameFile_msh)
%     disp('Done')
% end


disp('open GID FILE:')
%clipboard('copy',NAMEFILEOPEN)
disp(NAMEFILEOPEN)

DATA = DefaultField(DATA,'OPEN_GID',0) ;
if DATA.OPEN_GID ==1
    TTTT = ['gidpost ',NAMEFILEOPEN] ;
    unix(TTTT);
end