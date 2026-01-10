function [stressGLO,setElementsRED] = ...
    StressReconstruction(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim)

setElementsRED = [] ;
stressGLO = [] ;
%dbstop('246')
nDEF = size(BdomRED,2) ; 
%
if DATAINM.PRINT_STRESSES ==1  %& DATAINM.CUBATURE.ACTIVE ==1
    % Number of reduced elements per domain
    % ------------------------------------
    ngaus = size(posgp,2) ;
    if DATAINM.CUBATURE.ACTIVE ==1
        setElements =  large2small(setPoints,ngaus) ; % Selected elements
        %setGauss_red = small2large(setElements, ngaus) ;
    else
        setElements = [] ;
        
    end
   % DATA.setElementsRED = setElements ;
    %  rW = Wdom ;
    nstrain = size(BdomRED,1)/length(Wdom) ;
    rW = repmat(Wdom',nstrain,1) ;
    rW = rW(:) ;
    CBred = Cglo*BdomRED ;  % Celas*Bred
    nCOMPe = size(BdomRED,1) ;
    stressGLOv = zeros(nCOMPe*nDOM,1) ;
    nelem = length(Wdom)/ngaus ; 
    setElementsRED = [] ;
    for idom = 1:nDOM
        IND = (idom-1)*nDEF +1: idom*nDEF ;
        strainDOM_e = BdomRED*qDEF(IND) ;
        stressDOM_e =(Cglo*strainDOM_e)./rW ;
        IND = (idom-1)*nCOMPe+1:idom*nCOMPe ;
        stressGLOv(IND) = stressDOM_e ;
        if ~isempty(setElements)
            setLOC = (idom-1)*nelem + setElements ; 
            setElementsRED= [setElementsRED; setLOC] ; 
        end
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
    
    
end
