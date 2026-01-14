function  [stressGLO_glo,stressVONMISES_glo,MAXstressVONMISES] = ...
    FullOrderStresses(ientity,DATAIN,DATA_REFMESH,stressGLO_glo,stressVONMISES_glo,qDEF,DATAROM,ELEMS,nstrain,...
    SELECTED_DOMAINS,DATARUN,DATAadd)



% ---------------------
%% FULL-ORDER MODEL
%% --------------------
disp('Retrieving OFFLINE data')
tic
load(DATAIN.NAME_WS_MODES{ientity},'CgloDOM','Wdom','Bdom')
toc
disp('DONE')


if ~exist('DATARUN','var')
    DATARUN = [] ;
end

DATARUN = DefaultField(DATARUN,'USE_RECONSTRUCTION_MATRIX',0) ;
DATARUN = DefaultField(DATARUN,'ISNONLINEAR',0) ;

if ~isempty(DATAadd)
    stressST = DATAadd.stressST ;
    nDOMM = size(qDEF,2) ;
    stressST = reshape(stressST,[],nDOMM) ;
else
    stressST = [] ;
end


if isempty(stressST)
    if DATARUN.USE_RECONSTRUCTION_MATRIX == 0
        Bdom = Bdom*DATAROM{ientity}.BasisUdef ;
        CgloDOM = CgloDOM*Bdom ;
    else
        % C*B  = ReconstMat*(C*B)_{setPoints}
        CgloDOM = DATAROM{ientity}.HROMVAR.ReconsStresses*DATAROM{ientity}.HROMVAR.Celas_Bdom ;
    end
    
end




if DATAIN.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 0
    error(['Option not implemented'])
end
ngaus = size(DATA_REFMESH{ientity}.posgp,2) ;
Wst = repmat(Wdom',nstrain,1) ;
Wst = Wst(:) ;
nelemALL =  length(Wdom)/ngaus ;

DATAIN = DefaultField(DATAIN,'ONLY_SHOW_REDUCED_ELEMENTS',0) ; 
DATAIN = DefaultField(DATAIN,'ONLY_SHOW_REDUCED_ELEMENTS',0) ; 

setElementsShowREDUCED = DATAROM{ientity}.HROMVAR.setElements ; 
if DATAIN.ONLY_SHOW_REDUCED_ELEMENTS == 1 
    setElementsShow = setElementsShowREDUCED ; 
    nelem  = length(setElementsShow) ;
else
    setElementsShow = [] ; 
    nelem = nelemALL ; 
end  



if ~isempty(SELECTED_DOMAINS)
    ELEMS = ELEMS(SELECTED_DOMAINS{ientity});
end


stressGLO_glo{ientity} = zeros(nelem*nstrain,size(length(ELEMS),2)) ;
stressVONMISES_glo{ientity}= zeros(nelem,size(length(ELEMS),2)) ;
MAXstressVONMISES = zeros(1,size(length(ELEMS),2))  ;



if DATAIN.DO_NOT_COMPUTE_STRESSES == 0
    tic
    disp('Computing average stress on each FE element...')
    
    for ielem = 1:length(ELEMS)
        idom = ELEMS(ielem) ;
        disp(['Domain =',num2str(ielem)])
        if isempty(stressST)
        stressDOM_e = CgloDOM*qDEF(:,idom) ;  %
        else
            stressST_red = stressST(:,idom) ; 
            stressDOM_e = DATAROM{ientity}.HROMVAR.ReconsStresses*stressST_red ; 
        end
        
        if DATARUN.ISNONLINEAR == 0
            for istrain =1:nstrain
                stressDOM_e(istrain:nstrain:end) = stressDOM_e(istrain:nstrain:end)./Wdom ;   % Cglo already includes WEIGHTS
            end
        end
        [stressDOM_e] = AverageStressOnElements(stressDOM_e,Wst,nelemALL,nstrain,ngaus) ;
        
        % Von Mises
        stressDOM_ELEM = reshape(stressDOM_e,nstrain,[]) ;
        
        if ~isempty(setElementsShow)
             stressDOM_ELEM = stressDOM_ELEM(:,setElementsShow) ; 
        end
        
        [ stressVONMISES_e ] =  VonMises_Stress(stressDOM_ELEM) ;
        MAXstressVONMISES(ielem) = max(stressVONMISES_e) ;
        stressGLO_glo{ientity}(:,ielem) = stressDOM_ELEM(:) ;
        stressVONMISES_glo{ientity}(:,ielem) = stressVONMISES_e' ;
    end
    disp('...Done')
    toc
else
    stressVONMISES=[] ;
    stressGLO = [] ;
    MAXstressVONMISES = [] ;
    
end

end