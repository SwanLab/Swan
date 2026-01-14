function  [MAXstressVONMISES] = ...
    ReducedOrderStresses(ientity,DATAIN,DATA_REFMESH,stressGLO_glo,stressVONMISES_glo,...
    qDEF,DATAROM,ELEMS,nstrain,DATAadd)


if nargin == 0
    load('tmp2.mat')
end

% ---------------------
%% FULL-ORDER MODEL
%% --------------------
% disp('Retrieving OFFLINE data')
% tic
% load(DATAIN.NAME_WS_MODES{ientity},'CgloDOM','Wdom','Bdom')
% toc
% disp('DONE')
%%%%%%%%%%%%%%%%%%%%%%%%
%if ~isfield(DATAROM{ientity},'HROMVAR')
Celas_Bdom = DATAROM{ientity}.HROMVAR.Celas_Bdom;   % Product Celas times B matrix (reduced points )
%end

if DATAIN.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 0
    error(['Option not implemented'])
end
% ngaus = size(DATA_REFMESH{ientity}.posgp,2) ;
% Wst = repmat(Wdom',nstrain,1) ;
% Wst = Wst(:) ;
% nelem =  length(Wdom)/ngaus ;



% stressGLO_glo{ientity} = [] ;
% stressVONMISES_glo{ientity}= [] ;
MAXstressVONMISES = zeros(1,size(length(ELEMS),2))  ;

DATAIN = DefaultField(DATAIN,'POST_PROCESS_COARSE_STRESS_AVERAGE',0) ;  % Post-process coarse-scale average stresses
DATAIN = DefaultField(DATAIN,'POST_PROCESS_VON_MISES_PER_MATERIAL',0) ; % Set it to the number of material you wish to plot
                                                 % If 0, then it is the
                                                 % maximum, with no
                                                 % distinction of the
                                                 % material it comes from
                                                 % (18-Feb-2020)
                                                 
ientity = 1; 
MaterialType = DATA_REFMESH{ientity}.MaterialType ;
SetElemtReduced = DATAROM{ientity}.HROMVAR.setElements ; 
ElemSearchMax = 1:length(SetElemtReduced)  ; 


if DATAIN.POST_PROCESS_VON_MISES_PER_MATERIAL >0
    MaterialType = MaterialType(SetElemtReduced) ; 
   ElemSearchMax = find(MaterialType==DATAIN.POST_PROCESS_VON_MISES_PER_MATERIAL)    ; 
  
end

if DATAIN.DO_NOT_COMPUTE_STRESSES == 0
    tic
    disp('Computing average stress on each FE element')
    
    if ~isempty(DATAadd)
        stressST = DATAadd.stressST ;
        nDOMM = size(qDEF,2) ;
        stressST = reshape(stressST,[],nDOMM) ;
        
       COMPUTE_STRESS_AVERAGE_VIA_HOMOGENIZATION_OPERATOR = 1;   % iNTRO. 19-jAN-2020
            
          W = DATAROM{ientity}.HROMVAR.WdomRED ; % Weights
            VOL = sum(W) ;
            
       if COMPUTE_STRESS_AVERAGE_VIA_HOMOGENIZATION_OPERATOR == 0
            W = DATAROM{ientity}.HROMVAR.WdomRED ; % Weights
            VOL = sum(W) ;
           nstrain =  DATAROM{ientity}.HROMVAR.nstrain ;
            stressAVERAGE = zeros(nstrain,size(stressST,2)) ;
            for istrain = 1:nstrain
                ROWS = istrain:nstrain:size(stressST,1) ;
                stressSTloc = stressST(ROWS,:) ;
                wSTRESS = sum(bsxfun(@times,stressSTloc,W),1) ;
                stressAVERAGE(istrain,:) = wSTRESS/VOL ;
            end
       else
           stressAVERAGE  =  DATAROM{ientity}.HROMVAR.HOMOGENIZATION_OPERATOR*stressST/VOL ;
       end
            
       
        
        
        
    else
        stressST = [] ;
        stressAVERAGE = [] ;
    end
    
    
    if DATAIN.POST_PROCESS_COARSE_STRESS_AVERAGE == 0 || isempty(stressAVERAGE)
        for ielem = 1:length(ELEMS)
            idom = ELEMS(ielem) ;
            disp(['Domain =',num2str(ielem)])
            
            % Stresses at the selected Gauss points
            if isempty(stressST)
                stressDOM_e = Celas_Bdom*qDEF(:,idom) ;
            else
                stressDOM_e = stressST(:,idom) ;
            end
            stressDOM_ELEM = reshape(stressDOM_e,nstrain,[]) ;
            
            
            %         for istrain =1:nstrain
            %             stressDOM_e(istrain:nstrain:end) = stressDOM_e(istrain:nstrain:end)./Wdom ;   % Cglo already includes WEIGHTS
            %         end
            %     [stressDOM_e] = AverageStressOnElements(stressDOM_e,Wst,nelem,nstrain,ngaus) ;
            
            % Von Mises
            
            [ stressVONMISES_e ] =  VonMises_Stress(stressDOM_ELEM) ;
            [MAXstressVONMISES(ielem) INNDDDDD]= max(stressVONMISES_e(ElemSearchMax)) ;
            % stressGLO_glo{ientity}(:,ielem) = stressDOM_e ;
            % stressVONMISES_glo{ientity}(:,ielem) = stressVONMISES_e' ;
        end
        disp('...Done')
        toc
        
    else
      %  MAXstressVONMISES = stressAVERAGE ; 
         [ MAXstressVONMISES ] =  VonMises_Stress(stressAVERAGE(ElemSearchMax)) ;
        
    end
    
    
    
else
    %     stressVONMISES=[] ;
    %     stressGLO = [] ;
    MAXstressVONMISES = [] ;
    
end



end