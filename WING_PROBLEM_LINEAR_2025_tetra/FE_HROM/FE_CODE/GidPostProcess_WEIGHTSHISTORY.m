function [NameFile_res ]= GidPostProcess_WEIGHTSHISTORY(COOR,CN,TypeElement, ...
    posgp,NameFile_res,MaterialType,IND_ELEM_MESHES,...
    DATA_REFMESH,DATA,DATAIN);
% Post-processing of results using GID
%dbstop('5')
if nargin==0
    load('tmp1.mat')
end


% Writing results file
% -----------------------------------------------------
% Nodal vectors
% -----------------------------------------------------
NODAL_VECTOR = [] ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------
% STRESS DATA. 1D mesh
% -----------------------



nelem_TOT =  0 ;

imesh = 1 ;
MESH(imesh).TypeElement = TypeElement{imesh} ;
MESH(imesh).posgp = posgp{imesh} ;
MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
nelemE(imesh) = size(CN{imesh},2) ;
nelem_TOT = nelem_TOT +  size(CN{imesh},1) ;

%if ~isempty(STRESSDATA.MAX_VONMISES)
% scalar defined on GAUSS POINTS, mesh 1

MESH(imesh).GAUSS_SCALAR(1).NAME = 'WEIGHTS(perc.VOL)' ;
MESH(imesh).GAUSS_SCALAR(1).VAR =DATA.WEIGHTS ;
MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = DATA.ELEMENTS_WITH_WEIGHTS ;


ndim = size(COOR,2) ;
DATALOC = [] ;
DATAIN = DefaultField(DATAIN,'TIMELOC',1) ;
DATALOC.TIMELOC = DATAIN.TIMELOC;
NODAL_SCALAR = [] ;
OPERfe = [] ;



fid_res = fopen(NameFile_res,'wt');
fprintf(fid_res,'GiD Post Results File 1.0 \n');
nmesh = length(MESH) ;

for imesh = 1:length(MESH) ;
    MESHLOC = MESH(imesh) ;
    MESHLOC = DefaultField(MESHLOC,'posgp',[]) ;
    MESHLOC = DefaultField(MESHLOC,'TypeElement','Linear') ;
    
    if isempty(MESHLOC.posgp)
        npg = 1;
    else
        npg = size(MESHLOC.posgp,2);
    end
    fprintf(fid_res,['GaussPoints "',MESHLOC.NAMEMESH,'" Elemtype ',MESHLOC.TypeElement,'\n']);
    fprintf(fid_res,['Number of Gauss Points: ',num2str(npg),'\n']);
    fprintf(fid_res,'Nodes not included\n');
    fprintf(fid_res,'Natural Coordinates: Internal\n');
    fprintf(fid_res,'End GaussPoints\n');
    npG(imesh) = npg ;
end

  
DATA.TIME_DISCRETIZATION = 1:length(DATA.WEIGHTS) ; %  = DefaultField(DATA,'TIME_DISCRETIZATION',1) ;
DATA = DefaultField(DATA,'STEPS_TO_STORE',1:length(DATA.TIME_DISCRETIZATION)) ;
DATA = DefaultField(DATA,'STEPS_TO_PRINT', 1:length(DATA.STEPS_TO_STORE)) ;
DATA.STEPS_TO_PRINT = DATA.TIME_DISCRETIZATION ; 
 

for istepLOC = 1:length(DATA.STEPS_TO_PRINT)       
    istep = DATA.STEPS_TO_PRINT(istepLOC) ;
    time_step = DATA.TIME_DISCRETIZATION(DATA.STEPS_TO_STORE(istep)) ;
    
        
    % GAUSS RESULTS
    % --------------
    for imesh = 1:length(MESH)
        MESHloc =  MESH(imesh) ;  ;
        NAMEMESH_loc = MESHloc.NAMEMESH;
        npg = npG(imesh);
        MESHloc= DefaultField(MESHloc,'GAUSS_SCALAR',[]) ;
        NAMEFIELDS = [] ;
        for ivar = 1:length(MESHloc.GAUSS_SCALAR)
            NAMEFIELDS.VARIABLE = MESHloc.GAUSS_SCALAR(ivar).NAME ;
            VAR = MESHloc.GAUSS_SCALAR(ivar).VAR{istep}(:) ;
            VAR = VAR/sum(VAR)*100 ; 
            ELEMENTS = MESHloc.GAUSS_SCALAR(ivar).ELEMENTS{istep}' ;
            GidPostPlotScalarGaussNEW(VAR,npg,ELEMENTS,time_step,NAMEFIELDS,fid_res,NAMEMESH_loc) ;
        end
        
        
      
        
    end
    
    
    
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% STRESSES ........................
    %     % ...................................
    %     if ~isempty(stressGLO)
    %         NAMEFIELDS.VARIABLE = 'STRESS' ;
    %         NAMEFIELDS.COMP  = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'} ;
    %         setElementsRED = [] ;
    %         GidPostPlotMatGauss(stressGLO,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% REDUCED STRESSES ........................
    %     % ...................................
    %     DATA  =DefaultField(DATA,'setElementsRED',[]) ;
    %     if ~isempty(stressGLO) & ~isempty(DATA.setElementsRED);
    %         NAMEFIELDS.VARIABLE = 'GAPPY STRESS' ;
    %         NAMEFIELDS.COMP  = {'StressG-xx','StressG-yy','StressG-zz','StressG-xy','StressG-yz','StressG-xz'} ;
    %         setElementsRED = DATA.setElementsRED ;
    %         var = stressGLO(:,DATA.setElementsRED) ;
    %         GidPostPlotMatGauss(var,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% STRAINS ........................
    %     % ...................................
    %     if ~isempty(strainGLO)
    %         NAMEFIELDS.VARIABLE = 'STRAIN' ;
    %         NAMEFIELDS.COMP  = {'Strain-xx','Strain-yy','Strain-zz','Strain-xy','Strain-yz','Strain-xz'} ;
    %         setElementsRED = [] ;
    %         GidPostPlotMatGauss(strainGLO,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% VON MISES STRESS ........................
    %     % ...................................
    %     DATA = DefaultField(DATA,'PRINT_VONMISES_STRESS',1) ;
    %     DATA = DefaultField(DATA,'stressVONMISES',[]) ;
    %     if ~isempty(DATA.stressVONMISES) & DATA.PRINT_VONMISES_STRESS == 1
    %         NAMEFIELDS.VARIABLE = 'VON_MISES' ;
    %         setElementsRED = [] ;
    %         GidPostPlotScalarGauss(DATA.stressVONMISES,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %     DATA = DefaultField(DATA,'MAXstressVONMISES',[]) ;
    %     if ~isempty(DATA.MAXstressVONMISES)
    %         NAMEFIELDS.VARIABLE = 'MAX_VON_MISES' ;
    %         setElementsRED = [] ;
    %         GidPostPlotScalarGauss(DATA.MAXstressVONMISES,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %     if ~isempty(DATA.Forces1D)
    %
    %         fff = fieldnames(DATA.Forces1D) ;
    %         for ifield = 1:length(fff)
    %             FORCE = DATA.Forces1D.(fff{ifield}) ;
    %             if ~isempty(FORCE)
    %                 NAMEFIELDS.VARIABLE = fff{ifield} ;
    %                 GidPostPlotScalarNodes(FORCE,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
    %             end
    %         end
    %     end
    
end




fclose(fid_res);

