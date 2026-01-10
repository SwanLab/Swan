function [stressGLO,setElementsRED,stressVONMISES,MAXstressVONMISES] = ...
    StressReconstructionMULTI_elem(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim,...
    COL,COLloc)

if nargin == 0
    load('tmp.mat')
end
strainGLO = [] ;
stressGLO = [] ;
stressVONMISES = [] ;
MAXstressVONMISES = [] ;
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
        %setGauss_red = small2large(setElements, ngaus) ;
    else
        setElements = [] ;
        
    end
    % DATA.setElementsRED = setElements ;
    %  rW = Wdom ;
    nstrain = size(BdomRED{1},1)/length(Wdom) ;
    rW = repmat(Wdom',nstrain,1) ;
    rW = rW(:) ;
    %  for itype = 1:length(BdomRED)
    %     CBred = Cglo*BdomRED{itype} ;  % Celas*Bred
    [ nCOMPe nDEF]= cellfun(@size,BdomRED) ;
    %  end
    nCOMPe = nCOMPe(1) ;
    if DATAINM.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 0
        stressGLOv = zeros(nCOMPe*nDOM,1) ;
        
    else
        nCOMPe = nCOMPe/ngaus ;
        stressGLOv = zeros(nCOMPe*nDOM,1) ;
        stressVONMISES = zeros(nCOMPe/nstrain,1) ;
        MAXstressVONMISES = zeros(nDOM,1) ;
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
        
        
        
        IND = (idom-1)*nCOMPe+1:idom*nCOMPe ;
        stressGLOv(IND) = stressDOM_e ;
        if ~isempty(stressVONMISES)
            IND = (idom-1)*nCOMPe/nstrain+1:idom*nCOMPe/nstrain ;
            stressVONMISES(IND) = stressVONMISES_e ;
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
