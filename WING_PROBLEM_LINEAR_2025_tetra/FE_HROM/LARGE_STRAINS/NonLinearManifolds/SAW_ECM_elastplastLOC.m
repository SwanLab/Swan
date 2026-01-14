function [ECMdata_cluster,setCandidates,LIST_OF_CANDIDATES] = SAW_ECM_elastplastLOC(SNAPfint_lin,SNAPfint_non,DATA,wSTs,DATAoffline,qPLAST)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
    %     DATAoffline.ECM_Ratio_Number_Clusters_Snapshots = 1 ;
    DATAoffline.ECM_Number_Snapshots_Overlapping =1 ;
    DATAoffline.IncludeElastic_SAW_ECM =1;
       DATAoffline.IncludeConstantFunction_SAW_ECM = 1; 
        DATAoffline.errorFINT  =1e-4 ; 
       % DATAoffline.SetCandidates_given_by_the_userECM = [ 12440       13234       14769       15239       15796       15797       16325       16885       16886       16889       17153       17488       17492  17714       18119] ; 

    close all
    %     DATAoffline.errorFINT = 1e-6 ;
end


DATAoffline = DefaultField(DATAoffline,'ECM_Number_Snapshots_Overlapping',1) ;

% INTRODUCING OVERLAPPING  FOR NONLINEAR SNAPSHOTS

SNAPfint_non_OVERLAPPED = SNAPfint_non ;
Wfe = diag(sparse(wSTs)) ;
DATAloccc.TOL =0;
NumberBasisFnon_overlap = zeros(size(SNAPfint_non)) ;
for icluster= 1:length(SNAPfint_non)
    icluster_back = max(1,icluster-DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluster_forw = min(length(SNAPfint_non),icluster+DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluser_select = [icluster_back:icluster_forw] ;
    SNAPloc = cell2mat(SNAPfint_non(icluser_select)) ;
    
    if icluster == 1
        [SNAPfint_non_OVERLAPPED{icluster},SSS,~,sqW] = WSVDT(SNAPloc,Wfe,DATAloccc) ;
        DATAloccc.Mchol = sqW ;
    else
        [SNAPfint_non_OVERLAPPED{icluster},SSS,~,sqW] = WSVDT(SNAPloc,[],DATAloccc) ;
    end
    NumberBasisFnon_overlap(icluster) = length(SSS) ;
    
end

figure(1434)
hold on
xlabel('Number of time step (plastic loading)')
ylabel('Number of basis functions, nonlinear')
title(['Number of nonlinear basis functions, OVERLAPPING level =',num2str(DATAoffline.ECM_Number_Snapshots_Overlapping)]) ;
bar(NumberBasisFnon_overlap);


SNAPfint_non = SNAPfint_non_OVERLAPPED ;

% Basis matrix for internal forces
% **********************************
%if iscell(SNAPfint)
% Elastic modes ---how many elastic modes are here?
DATAoffline = DefaultField(DATAoffline,'IncludeElastic_SAW_ECM',1)  ;

if  DATAoffline.IncludeElastic_SAW_ECM  ==1
    SNAPfint_lin = cell2mat(SNAPfint_lin) ;
    TOL = 1e-6;
    [PhiELAS,S,V] = SRSVD(SNAPfint_lin,TOL) ;
    %Wfe = diag(sparse(wSTs)) ;
    [PhiELAS] = WSVDT(PhiELAS,[],DATAloccc) ;
    disp(['Number of elastic modes (internal forces) = ',num2str(size(PhiELAS,2))]);
    %  Enlarged to include the constant function
    %DATAloccc.Mchol = sqW ;
    
    DATAoffline = DefaultField(DATAoffline,'IncludeConstantFunction_SAW_ECM',1)  ;
    
    if DATAoffline.IncludeConstantFunction_SAW_ECM == 1
        PhiELAS = WSVDT([PhiELAS,ones(size(PhiELAS,1),1) ],[],DATAloccc ) ;
    end
    %%%% Next we compute the orthogonal complement
    SNAPfint_w = SNAPfint_non;
    TOL = DATAoffline.errorFINT ;
    NumberModesCluster = zeros(size(SNAPfint_non)) ;
    for isnap = 1:length(SNAPfint_non)
        [UU,SS,~] =   SRSVD({sqW*SNAPfint_non{isnap},sqW*PhiELAS},[TOL,0]);
        SNAPfint_w{isnap} =  UU ;
        NumberModesCluster(isnap) = length(SS) ;
    end
    
    figure(42)
    hold on
    xlabel('Snapshot')
    ylabel('Number of internal force modes')
    bar(NumberModesCluster)
    
else
     TOL = DATAoffline.errorFINT ;
       NumberModesCluster = zeros(size(SNAPfint_non)) ;
      for isnap = 1:length(SNAPfint_non)
        [UU,SS,~] =   SRSVD(sqW*SNAPfint_non{isnap},TOL);
        SNAPfint_w{isnap} =  UU ;
        NumberModesCluster(isnap) = length(SS) ;
    end
     
    figure(42)
    hold on
    xlabel('Snapshot')
    ylabel('Number of internal force modes')
    bar(NumberModesCluster)
end

%
[ECMdata_cluster,setCandidates,LIST_OF_CANDIDATES] = ...
    SAW_ECM(wSTs,DATAoffline,SNAPfint_w,DATA) ;

disp(['Number of Candidate Points'])
disp(length(setCandidates))
disp(['setCandidates = '])
disp(setCandidates')



% % %
% % %
% % %
% % %
% % %
% % % DATAoffline = DefaultField(DATAoffline,'ListElementsExclude_fromGID',[]) ;
% % % if ~isempty(DATAoffline.ListElementsExclude_fromGID)
% % %     ListElementsToExclude = load(DATAoffline.ListElementsExclude_fromGID) ;
% % %     ListElementsToExclude = ListElementsToExclude(:,1) ;
% % %     ngausELEM = DATA.MESH.ngaus_STRESS;
% % %     ListGaussToExclude = small2large(ListElementsToExclude,ngausELEM) ;
% % %     INDSEL  = setdiff(1:size(Q,1),ListGaussToExclude) ;
% % % else
% % %     INDSEL  = 1:size(Q,1) ;
% % % end
% % %
% % % DATA_ECM.IND_POINTS_CANDIDATES = INDSEL ;
% % % DATAoffline = DefaultField(DATAoffline,'errorECM',0) ;
% % % DATA_ECM.TOL = DATAoffline.errorECM ;
% % % DATAoffline = DefaultField(DATAoffline,'USE_SELECTIVE_DECM',1) ; % = 1;
% % % if DATAoffline.USE_SELECTIVE_DECM == 0
% % %     % Version before 3-DEc-2021
% % %     [setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;
% % % else
% % %     % New version (after 3-DEc-2021)
% % %     [setPoints,wRED,ERROR_GLO,DATAOUT]= DiscreteEmpiricalCubatureMethod(Q',wSTs,DATA_ECM)  ;
% % %
% % % end
% % %
% % % DATAOUT.TIMEdecm = toc(TIMEdecm) ;
% % %
% % % disp(['Time to select the integration points =',num2str(DATAOUT.TIMEdecm)])
% % %
% % %
% % % SNAPfint = cell2mat(SNAPfint) ;
% % %
% % % IntegralExact =SNAPfint'*wSTs  ;
% % % IntegrationError = IntegralExact - SNAPfint(setPoints,:)'*wRED ;
% % % IntegrationError = norm(IntegrationError)/norm(IntegralExact);
% % % disp(['Actual integration error using DECM= ',num2str(IntegrationError*100),' % (prescribed tolerance fint =',num2str(DATAoffline.errorFINT*100), '%'])