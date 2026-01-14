function [ECMdata,BasisPone,Spone] = Stresses_EmpiricalCubMethod(NAMEsnap_base,DATAoffline,BasisU,CASES,DOFr,DOFl,DATAparamSTUDY,OPERFE,MATPRO,DATA)


% Now we have to compute the PK2 stresses produced by these strains, and check
% that all blocks are approximated to an accuracy theshold
% DATAoffline.\errorSTRESS
[BasisPone,Spone,nmodesS] = StressesBasisFromDisplacements(BasisU,CASES,NAMEsnap_base,DOFr,DOFl,DATAparamSTUDY,OPERFE,MATPRO,DATAoffline,DATA) ;
% -------------------
% HYPERREDUCTION
% -------------------
%
% else
BstRED_l = OPERFE.Bst(:,DOFl)*BasisU ;
%
% end
if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
    % ******************
    % Discrete ECM
    % ******************
    DATA_GENGAUSS.ACTIVE =  0 ; 
    if DATA_GENGAUSS.ACTIVE == 0
        [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM(BstRED_l,BasisPone,DATA,OPERFE.wSTs,DATAoffline) ;
        ECMdata.setPoints = setPoints ;
        ECMdata.wRED = wRED ;
        proporP = length(setPoints)/DATA.MESH.ngausT*100;
        disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])
        
        setElements = large2smallREP(setPoints,DATA.MESH.ngaus) ;
        disp('****************************+')
        disp(['List of selected m = ',num2str(length(setElements)),' elements'])
        disp(num2str(setElements'))
        %clipboard('copy',num2str(setElements'));
        ECMdata.setElements = setElements ;
        
    else
        % ****************************************************************
        % Function for computing the position of the integration points
        % ****************************************************************
        Nst = OTHER_output.Nst ;
        NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
        if  ~exist(NAME_MODES_FOLDER)
            mkdir(NAME_MODES_FOLDER)
        end
        DATA_GENGAUSS.NameFileMesh_ECM = [NAME_MODES_FOLDER,NAME_BASE,'DECMpoints'] ;
        DATALOC = [] ;
        
        DATA_GENGAUSS.NameFileMesh_FINT = [NAME_MODES_FOLDER,NAME_BASE,'InternalForceModes'] ;
        DATALOC = [] ;
        
        DATA_GENGAUSS.NameFileMesh_CECM = [NAME_MODES_FOLDER,NAME_BASE,'CECMpoints'] ;
        DATALOC = [] ;
        
        [ECMdata] = ContinuousECM(BstRED_l,BasisPone,DATA,OPERFE.wSTs,DATAoffline,DATA_GENGAUSS,...
            MESH,Nst) ;
    end
else
    SNAPredFINT = BasisF_from_BasisStress_PK1_ELEMS(BstRED_l,BasisPone,DATA, OPERFE.wSTs)  ;
    wSTs_LOC = ones(size(SNAPredFINT,1),1) ;
    %     sqrt_wST = sqrt(OPERFE.wSTs) ;
    %     SNAPredFINT = bsxfun(@times,SNAPredFINT,sqrt_wST) ;
    % Determine an  orthogonal basis matrix $Q$ for the column space of $\SNAPredFINTw{}{}$
    %DATAoffline.errorFINT = 1e-3;
    DATAsvd.RELATIVE_SVD = 1;
    [Q,S,V,eSVD,Rsup] = RSVDT(SNAPredFINT,DATAoffline.errorFINT,[],0,DATAsvd) ;
    if DATAoffline.errorFINT == 0
        ifig = 3000 ;
        SVDplotERROR_local(S,ifig) ;
    end
    % % Enlarge the basis matris for SNAPredFINT
    a  = wSTs_LOC - Q*(Q'*wSTs_LOC) ;
    if norm(a) > 1e-10
        a = a/norm(a) ;
        Q = [Q,a] ;
    end
    % Empirical cubature method
    % -------------------------
    DATA_ECM = [] ;
    DATA_ECM.TOL = DATAoffline.errorECM ;
    [setElements,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs_LOC,DATA_ECM)  ;
    disp(['Element-based approach *************************++'])
    disp(['List of selected m = ',num2str(length(setElements)),' elements'])
    disp(num2str(setElements'))
    figure(345)
    hold on
    xlabel('Points')
    ylabel('Weights')
    
    bar(sort(wRED,'descend'))
    
    % Determine set of points
    setPoints = small2large(setElements,DATA.MESH.ngaus_STRESS) ;
    disp(['Total number of Gauss points = ',num2str(length(setPoints))])
    wRED = repmat(wRED',DATA.MESH.ngaus_STRESS,1) ;
    wRED = wRED(:).*OPERFE.wSTs(setPoints,:) ;
    
    ECMdata.setPoints = setPoints ;
    ECMdata.wRED = wRED ;
    ECMdata.setElements = setElements ;
    
end

disp(['*********************************************************************'])
disp(['Number of displacement modes = ',num2str(size(BasisU,2))])
disp(['Number of PK2-stress modes = ',num2str(nmodesS)])
disp(['Number of PK1-stress modes = ',num2str(size(BasisPone,2))])