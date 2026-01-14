function [stressGLO_glo,STRESS_DATA,SELECTED_DOMAINS] = ...
    StressDom_SlicesJoints(DATAROM,MESH1D,DATAIN,qDEF,DATA_REFMESH,NO_STRESSES,DATARUN,DATAadd)

if nargin == 0
    load('tmp2.mat')
end
DATAIN = DefaultField(DATAIN,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',1) ;
DATAIN = DefaultField(DATAIN,'DO_NOT_COMPUTE_STRESSES',0) ;
if size(MESH1D.COOR,2) == 3
    nstrain = 6 ;
else
    nstrain = 4;
end
nentities = length(DATAROM) ;
stressGLO_glo = cell(nentities,1) ;
nelem1D = size(MESH1D.CN,1) ;
stressVONMISES_glo = cell(nentities,1) ;
STRESS_DATA.MAX_VONMISES = zeros(nelem1D,1) ;
DATAIN = DefaultField(DATAIN,'COMPUTE_STRESSES_AT_REDUCED_POINTS',1) ;
SELECTED_DOMAINS = cell(nentities,1) ;
for ientity = 1:length(DATAROM)
    switch  MESH1D.PROP(ientity).TYPE
        case 'JOINT'
            ISJOINT = 1;
        otherwise
            ISJOINT = 0 ;
    end
    ELEMS = find(MESH1D.MaterialType == ientity) ;
    
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
            SelectingDomains(DATAIN,ISJOINT,ientity,MAXstressVONMISES,NO_STRESSES,...
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
            ELEMS,nstrain,[]) ;
        % ----------------------------------
        
        % SELECTING THE DOMAINS WHOSE STRESS FIELD WILL BE SHOWN
        % ---------------------------------------------------------------
        [SELECTED_DOMAINS] =  ...
            SelectingDomains(DATAIN,ISJOINT,ientity,MAXstressVONMISES,NO_STRESSES,SELECTED_DOMAINS,ELEMS) ;
        stressGLO_glo{ientity} =  stressGLO_glo{ientity}(:,SELECTED_DOMAINS{ientity}) ;
        stressVONMISES_glo{ientity} =  stressVONMISES_glo{ientity}(:,SELECTED_DOMAINS{ientity}) ;
    end
    
    
    %%% GID PRINTING
    if ~isempty(stressGLO_glo{ientity})
        if size(MESH1D.COOR,2) == 3
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

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%



function [SELECTED_DOMAINS] =  ...
    SelectingDomains(DATAIN,ISJOINT,ientity,MAXstressVONMISES,NO_STRESSES,...
    SELECTED_DOMAINS,ELEMS)

nmax = length(ELEMS);

if ~isempty(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER) & DATAIN.DO_NOT_COMPUTE_STRESSES ==0  & ISJOINT ==0
    % error('Modify this part of the code  !!!!!')
    NUMBERprint = min(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER,nmax) ;
    switch DATAIN.DOMAINS_POSTPROCESS_SELECT.VARIABLE
        case 'VONMISES'
            [~,INDEXES] =  sort(MAXstressVONMISES,'descend') ;
            NUMBERprint = min(NUMBERprint,length(INDEXES)) ;
            SELECTED_DOMAINS_loc = unique(INDEXES(1:NUMBERprint)) ;
            SELECTED_DOMAINS_loc = [1,SELECTED_DOMAINS_loc,length(MAXstressVONMISES)] ; % First and last domains are included
            SELECTED_DOMAINS_loc = unique(SELECTED_DOMAINS_loc) ;
            SELECTED_DOMAINS{ientity} = SELECTED_DOMAINS_loc ;
            
            %             if  ~isempty( stressVONMISES_glo{ientity} )
            %             stressVONMISES_glo{ientity}  =  stressVONMISES_glo{ientity}(:,SELECTED_DOMAINS_loc) ;
            %             stressGLO_glo{ientity} =  stressGLO_glo{ientity}(:,SELECTED_DOMAINS_loc) ;
            %             end
            %      DISP3D = DISP3D(:,SELECTED_DOMAINS) ;
            %    DATAIN.DOMAINS_POSTPROCESS = SELECTED_DOMAINS ;
        otherwise
            error('Option not implemented')
    end
else
    if ~isempty(DATAIN.DOMAINS_POSTPROCESS) && NO_STRESSES == 0 && DATAIN.DO_NOT_COMPUTE_STRESSES ==0
        
        [DDD iA iB]= intersect(DATAIN.DOMAINS_POSTPROCESS,ELEMS) ;
        
        % stressVONMISES_glo{ientity}  =  stressVONMISES_glo{ientity}(:,iB) ;
        % stressGLO_glo{ientity}  =  stressGLO_glo{ientity}(:,iB) ;
        SELECTED_DOMAINS{ientity} = iB ;
    elseif  isempty(DATAIN.DOMAINS_POSTPROCESS) && NO_STRESSES == 1
        % stressVONMISES =  [] ;
        % stressGLO =  [] ;
        SELECTED_DOMAINS = [1:length(ELEMS)];
    end
    
    
end

end



