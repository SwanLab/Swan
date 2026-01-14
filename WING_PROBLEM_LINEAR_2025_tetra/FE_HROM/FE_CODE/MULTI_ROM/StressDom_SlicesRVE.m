function [stressGLO_glo,STRESS_DATA,SELECTED_DOMAINS] = ...
    StressDom_SlicesRVE(DATAROM,MESH2D,DATAIN,qDEF,DATA_REFMESH,NO_STRESSES,DATAadd,DATARUN)

if nargin == 0
    load('tmp.mat')
end
DATAIN = DefaultField(DATAIN,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',1) ;
DATAIN = DefaultField(DATAIN,'DO_NOT_COMPUTE_STRESSES',0) ;
nMODESrb = size(DATAROM{1}.BasisIntRB{1},2) ;
nstrain = 6 ;

if nMODESrb == 3
    nstrain = 4 ;
end
nentities = length(DATAROM) ;
stressGLO_glo = cell(nentities,1) ;
nelem2D = size(MESH2D.CN,1) ;
stressVONMISES_glo = cell(nentities,1) ;
STRESS_DATA.MAX_VONMISES = zeros(nelem2D,1) ;
DATAIN = DefaultField(DATAIN,'COMPUTE_STRESSES_AT_REDUCED_POINTS',1) ;
SELECTED_DOMAINS = cell(nentities,1) ;
for ientity = 1:length(DATAROM)
    
    ELEMS = find(MESH2D.MaterialType == ientity) ;
    ISJOINT = 0 ;
    SELECTED_DOMAINS{ientity} = 1:length(ELEMS) ;
    
    if ~isfield(DATAROM{ientity},'HROMVAR')
        DATAIN.COMPUTE_STRESSES_AT_REDUCED_POINTS = 0;
    end
    
    
    if DATAIN.COMPUTE_STRESSES_AT_REDUCED_POINTS  == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computing REDUCED-order stresses (AT REDUCED-ORDER POINTS )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [MAXstressVONMISES] ...
            =  ReducedOrderStresses(ientity,DATAIN,DATA_REFMESH,stressGLO_glo,stressVONMISES_glo,qDEF,DATAROM,...
            ELEMS,nstrain,DATAadd) ;
        % SELECTING THE DOMAINS WHOSE STRESS FIELD WILL BE SHOWN
        % ---------------------------------------------------------------
        
        SELECTED_DOMAINS =  ...
            SelectingDomainsRVE(DATAIN,ISJOINT,ientity,MAXstressVONMISES,NO_STRESSES,...
            SELECTED_DOMAINS,ELEMS) ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computing full-order stresses (at all gauss points )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [stressGLO_glo,stressVONMISES_glo,MAXstressVONMISES_fullorder] ...
            =  FullOrderStresses(ientity,DATAIN,DATA_REFMESH,stressGLO_glo,stressVONMISES_glo,qDEF,DATAROM,...
            ELEMS,nstrain,SELECTED_DOMAINS,DATARUN,DATAadd) ;
        
        MAXVON_HYPER = max(MAXstressVONMISES) ;
        MAXVON_FE = max(MAXstressVONMISES_fullorder) ;
        DIFFF = norm(MAXVON_HYPER-MAXVON_FE)/norm(MAXVON_FE)*100 ;
        disp(['ERROR in predicting max. value VON MISES (%)=',num2str(norm(DIFFF))]) ;
        
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computing full-order stresses (at all gauss points )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [stressGLO_glo,stressVONMISES_glo,MAXstressVONMISES] ...
            =  FullOrderStresses(ientity,DATAIN,DATA_REFMESH,stressGLO_glo,stressVONMISES_glo,qDEF,DATAROM,...
            ELEMS,nstrain,[],DATARUN,DATAadd) ;
        % ----------------------------------
        
        % SELECTING THE DOMAINS WHOSE STRESS FIELD WILL BE SHOWN
        % ---------------------------------------------------------------
        [SELECTED_DOMAINS] =  ...
            SelectingDomainsRVE(DATAIN,ISJOINT,ientity,MAXstressVONMISES,NO_STRESSES,SELECTED_DOMAINS,ELEMS) ;
        stressGLO_glo{ientity} =  stressGLO_glo{ientity}(:,SELECTED_DOMAINS{ientity}) ;
        stressVONMISES_glo{ientity} =  stressVONMISES_glo{ientity}(:,SELECTED_DOMAINS{ientity}) ;
    end
    
    
    %%% GID PRINTING
    if ~isempty(stressGLO_glo{ientity})
        %  indGID = [1 2 3 6 4 5] ;
        if size(MESH2D.COOR,2) == 3
            indGID = [1 2 3 6 4 5] ;
        else
            indGID = [1 2 3 4] ;
        end
        
        stressGLO_old = stressGLO_glo{ientity};
        for iglo = 1:length(indGID)
            stressGLO_glo{ientity}(iglo:nstrain:end)  =  stressGLO_old(indGID(iglo):nstrain:end) ;
        end
        ngaus = 1;
        stressGLO_glo{ientity} =  reshape(stressGLO_glo{ientity},ngaus*nstrain,[]) ;
        STRESS_DATA.VONMISES{ientity} =  reshape(stressVONMISES_glo{ientity},ngaus,[]) ;
        STRESS_DATA.MAX_VONMISES(ELEMS) =   MAXstressVONMISES' ;
        
    else
        STRESS_DATA.VONMISES{ientity}  = [] ;
        STRESS_DATA.MAX_VONMISES  = [] ;
    end
    
    
end


end