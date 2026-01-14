function [ECMdata,setCandidates ]= ECM_clustersNCns(wSTs,DATAoffline,Q,DATA,setCandidates)

 
%
% end

DATAoffline = DefaultField(DATAoffline,'TOL_ECM',0) ; 

DATAoffline = DefaultField(DATAoffline,'NITERATIONS_NO_MORE_POINTS_negative_iterationECM',10) ; 
DATA.NITERATIONS_NO_MORE_POINTS_negative_iterationECM = DATAoffline.NITERATIONS_NO_MORE_POINTS_negative_iterationECM ; 
%   DATAoffline.
% if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
    % ******************
    % Discrete ECM
    % ******************
    % Q = QbasisMatrixIntegrand(BstRED_l,BasisPone,DATA,wSTs,DATAoffline) ;
    
    % Empirical cubature method (using as candidate points setCandidates)
    % -------------------------
    DATA_ECM = [] ;
    DATA_ECM.TOL = DATAoffline.TOL_ECM;
     if isempty(setCandidates)
        DATA_ECM.IND_POINTS_CANDIDATES = 1:length(wSTs) ; 
    else
        DATA_ECM.IND_POINTS_CANDIDATES = setCandidates ; 
    end
    
    [setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_CANDcompl(Q,[],wSTs,DATA_ECM)  ;
    
    setCandidates = unique([setCandidates;setPoints]) ; 
    
    
    %  [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM(BstRED_l,BasisPone,DATA,OPERFE.wSTs,DATAoffline) ;
    ECMdata.setPoints = setPoints ;
    ECMdata.wRED = wRED ;
    if any(isnan(wRED))
        save('ErrorW.mat','Q','wSTs','DATA_ECM')
        error('NaN weight !, see ErrorW.mat  ')
    end
    proporP = length(setPoints)/DATA.MESH.ngausT*100;
    disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])
    setElements = large2smallREP(setPoints,DATA.MESH.ngaus) ;
    disp('****************************+')
    disp(['List of selected m = ',num2str(length(setElements)),' elements'])
    disp(num2str(setElements'))
   % clipboard('copy',num2str(setElements'));
    ECMdata.setElements = setElements ;
    disp(['Number of total points = ',num2str(length(setCandidates))]) ;
    
      setElements_cand = large2smallREP(setCandidates,DATA.MESH.ngaus) ;
    disp('****************************+')
    disp(['List of selected candidates elements = ',num2str(length(setElements_cand)),' elements'])
    disp(num2str(setElements_cand'))
    
%     
% else
%     
%     error('Option not implemented yet')
%     SNAPredFINT = BasisF_from_BasisStress_PK1_ELEMS(BstRED_l,BasisPone,DATA, OPERFE.wSTs)  ;
%     wSTs_LOC = ones(size(SNAPredFINT,1),1) ;
%     %     sqrt_wST = sqrt(OPERFE.wSTs) ;
%     %     SNAPredFINT = bsxfun(@times,SNAPredFINT,sqrt_wST) ;
%     % Determine an  orthogonal basis matrix $Q$ for the column space of $\SNAPredFINTw{}{}$
%     %DATAoffline.errorFINT = 1e-3;
%     DATAsvd.RELATIVE_SVD = 1;
%     [Q,S,V,eSVD,Rsup] = RSVDT(SNAPredFINT,DATAoffline.errorFINT,[],0,DATAsvd) ;
%     if DATAoffline.errorFINT == 0
%         ifig = 3000 ;
%         SVDplotERROR_local(S,ifig) ;
%     end
%     % % Enlarge the basis matris for SNAPredFINT
%     a  = wSTs_LOC - Q*(Q'*wSTs_LOC) ;
%     if norm(a) > 1e-10
%         a = a/norm(a) ;
%         Q = [Q,a] ;
%     end
%     % Empirical cubature method
%     % -------------------------
%     DATA_ECM = [] ;
%     DATA_ECM.TOL = DATAoffline.errorECM ;
%     [setElements,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs_LOC,DATA_ECM)  ;
%     disp(['Element-based approach *************************++'])
%     disp(['List of selected m = ',num2str(length(setElements)),' elements'])
%     disp(num2str(setElements'))
%     figure(345)
%     hold on
%     xlabel('Points')
%     ylabel('Weights')
%     bar(sort(wRED,'descend'))
%     % Determine set of points
%     setPoints = small2large(setElements,DATA.MESH.ngaus_STRESS) ;
%     disp(['Total number of Gauss points = ',num2str(length(setPoints))])
%     wRED = repmat(wRED',DATA.MESH.ngaus_STRESS,1) ;
%     wRED = wRED(:).*OPERFE.wSTs(setPoints,:) ;
%     ECMdata.setPoints = setPoints ;
%     ECMdata.wRED = wRED ;
%     ECMdata.setElements = setElements ;
% end
% 
% if  isempty(ECMdata.setPoints)
%     BasisStwoZ =  InterpolationGaussVariablesECM(BasisStwo,ECMdata,DATA.MESH.ngaus_STRESS,DATA.MESH.nstrain) ;
% else
%     setIndices = small2large(ECMdata.setPoints,DATA.MESH.nstrain) ;
%     BasisStwoZ = BasisStwo(setIndices,:) ;
% end
% 
% ECMdata.coeff_reconstr_STWO =  (BasisStwoZ'*BasisStwoZ)\BasisStwoZ' ;




