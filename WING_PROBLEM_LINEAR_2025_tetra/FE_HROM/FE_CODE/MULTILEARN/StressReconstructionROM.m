function [stressGLO,setElementsRED,stressVONMISES,MAXstressVONMISES,DOMAINS_TO_INCLUDE] = ...
    StressReconstructionROM(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim,...
    COL,COLloc)

if nargin == 0
    load('tmp.mat')
end
stressGLO = cell(nDOM,1) ;
stressVONMISES = cell(nDOM,1) ;
MAXstressVONMISES = zeros(nDOM,1) ;
%dbstop('246')
%nDEF = size(BdomRED },2) ;
%


if DATAINM.PRINT_STRESSES ==1  %& DATAINM.CUBATURE.ACTIVE ==1
    % Number of reduced elements per domain
    % ------------------------------------
    ngaus = size(posgp,2) ;
    if DATAINM.CUBATURE.ACTIVE ==1
        setElements = cell(size(setPoints)) ;
        for itype = 1:length(setPoints)
            setElements{itype} =  large2small(setPoints{itype},ngaus) ; % Selected elements
            clipboard('copy',num2str(setElements{itype}'))
        end
    else
        setElements = [] ;
        
    end
    nstrain = size(BdomRED{1},1)/length(Wdom) ;
    rW = repmat(Wdom',nstrain,1) ;
    rW = rW(:) ;
    [ nCOMPe nDEF]= cellfun(@size,BdomRED) ;
    %  stressGLOv = cell(nDOM,1) ;
    if DATAINM.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 0
        nCOMPe = nCOMPe(1) ;
    else
        nCOMPe = nCOMPe/ngaus ;
    end
    %  wALL = zeros(nCOMPe*nDOM,1) ;
    nelem = length(Wdom)/ngaus ;
    setElementsRED = [] ;
    iacum = 0 ;
    
    for idom = 1:nDOM
        itype = COLloc(idom) ;
        IND = iacum +(1:nDEF(itype)) ; iacum = iacum + length(IND) ;
        strainDOM_e = BdomRED{itype}*qDEF(IND) ;
        stressDOM_e =(Cglo*strainDOM_e)./rW ;
        if DATAINM.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 1
            [stressDOM_e] = AverageStressOnElements(stressDOM_e,rW,nelem,nstrain,ngaus) ;
            % Von Mises
            stressDOM_ELEM = reshape(stressDOM_e,nstrain,[]) ;
            [ stressVONMISES_e ] =  VonMises_Stress(stressDOM_ELEM) ;
            MAXstressVONMISES(idom) = max(stressVONMISES_e) ;
        end
        
        
        
        %  IND = (idom-1)*nCOMPe+1:idom*nCOMPe ;
        stressGLO{idom} = stressDOM_e ;
        if ~isempty(stressVONMISES)
            % IND = (idom-1)*nCOMPe/nstrain+1:idom*nCOMPe/nstrain ;
            stressVONMISES{idom} = stressVONMISES_e' ;
        end
        
        %  wALL(IND) = rW ;
        if ~isempty(setElements)
            
            setLOC = (idom-1)*nelem + setElements{itype} ;
            setElementsRED= [setElementsRED; setLOC] ;
        end
        
        
        
    end
    if DATAINM.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 1
        ngaus = 1 ;
        
    end
    
    %%%%%%%%%%%%%%%%%
    if DATAINM.GID_1D_3D_print.ACTIVE == 1
        
        switch DATAINM.GID_1D_3D_print.METHOD_SELECTION
            case 'MAX_VON_MISES'
                ndom_to_include = DATAINM.GID_1D_3D_print.NUMBER_DOMAINS  ;
                [MAXV ia]= sort(MAXstressVONMISES,'descend');
                ndom_to_include = min(ndom_to_include,length(ia));
                DOMAINS_TO_INCLUDE = sort(ia(1:ndom_to_include)) ;
            case 'GIVEN'
                DOMAINS_TO_INCLUDE = DATAINM.GID_1D_3D_print.DOMAINS_TO_INCLUDE ;
        end
        
        stressGLOv = cell2mat(stressGLO(DOMAINS_TO_INCLUDE)) ;
        stressVONMISES = cell2mat(stressVONMISES(DOMAINS_TO_INCLUDE)) ;
        
    else
        stressGLOv = cell2mat(stressGLO ) ;
        stressVONMISES = cell2mat(stressVONMISES ) ;
        DOMAINS_TO_INCLUDE = [];
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%
    
    %%% GID PRINTING
    if ndim  == 3
        indGID = [1 2 3 6 4 5] ;
    else
        indGID = [1 2 3 ] ;
    end
    stressGLO = stressGLOv;
    for iglo = 1:length(indGID)
        stressGLO(iglo:nstrain:end)  =  stressGLOv(indGID(iglo):nstrain:end) ;
    end
    stressGLO =  reshape(stressGLO,ngaus*nstrain,[]) ;
    
else
    
    stressGLO= [] ;
    setElementsRED = [] ;
end

clipboard('copy',num2str(setElementsRED'))
