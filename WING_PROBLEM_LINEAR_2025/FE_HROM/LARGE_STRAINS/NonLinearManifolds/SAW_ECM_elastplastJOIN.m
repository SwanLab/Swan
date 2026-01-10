function  ECMdata = SAW_ECM_elastplastJOIN(DATAoffline,SNAPfint,DATA,OPERFE,qPLAST)

 
 
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp3.mat')
    DATAoffline.IntegrateVolumeExactlyECM = 0 ; 
end


DATAoffline = DefaultField(DATAoffline,'ECM_Number_Snapshots_Overlapping',1) ;
DATAoffline = DefaultField(DATAoffline,'IntegrateVolumeExactlyECM',1) ;

% Wfe = diag(sparse(OPERFE.wSTs)) ;
%  [PhiELAS,ssss,~,sqW] = SVDT(SNAPfint{1},Wfe) ;
%  disp(['Number of elastic modes (internal forces) = ',num2str(size(PhiELAS,2))]);
%    Enlarged to include the constant function
%  DATAloccc.Mchol = sqW ;
%
% DATAoffline = DefaultField(DATAoffline,'IncludeConstantFunction_SAW_ECM',1)  ;
%
% if DATAoffline.IncludeConstantFunction_SAW_ECM == 1
%PhiELAS = WSVDT([PhiELAS,ones(size(PhiELAS,1),1) ],[],DATAloccc ) ;
%end
%%%% Next we compute the orthogonal complement
SNAPfint_w = SNAPfint; 
DATAddd.ISRELATIVE = 1; 
sqW = sqrt(OPERFE.wSTs) ; 
sqWfe = diag(sparse(sqW)) ; 
NumberModesCluster = zeros(size(SNAPfint_w)) ; 
 
for icluster = 1:length(SNAPfint)
 
    disp(['icluster = ',num2str(icluster)])
    icluster_back = max(1,icluster-DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluster_forw = min(length(SNAPfint),icluster+DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluster_select = [icluster_back:icluster_forw] ;
     
        SNAPlocI = cell2mat(SNAPfint(icluster_select)) ;
        % INCLUDE ELASTIC SNAPSHOTS/integ. exact volum
        if DATAoffline.IntegrateVolumeExactlyECM == 1
        SNAPlocI = {[sqWfe*SNAPfint{1},sqrt(OPERFE.wSTs)],sqWfe*SNAPlocI} ;
        else
             SNAPlocI = {[sqWfe*SNAPfint{1}],sqWfe*SNAPlocI} ;
        end
        
        TOL = [0,DATAoffline.errorFINT] ; 
      DDD.HIDE_OUTPUT = 1; 
       [UU,SS,VV] = SRSVD(SNAPlocI,TOL,DDD) ; 
    
    SNAPfint_w{icluster} =  UU ;
    NumberModesCluster(icluster) = length(SS) ;
end

figure(42)
hold on
xlabel('Snapshot')
ylabel('Number of internal force modes')
bar(NumberModesCluster)
%
% else
%      TOL = DATAoffline.errorFINT ;
%        NumberModesCluster = zeros(size(SNAPfint_non)) ;
%       for isnap = 1:length(SNAPfint_non)
%         [UU,SS,~] =   SRSVD(sqW*SNAPfint_non{isnap},TOL);
%         SNAPfint_w{isnap} =  UU ;
%         NumberModesCluster(isnap) = length(SS) ;
%     end
%
%     figure(42)
%     hold on
%     xlabel('Snapshot')
%     ylabel('Number of internal force modes')
%     bar(NumberModesCluster)
% end

%
[ECMdata_cluster,setCandidates,LIST_OF_CANDIDATES] = ...
    SAW_ECM(OPERFE.wSTs,DATAoffline,SNAPfint_w,DATA) ;

disp(['Number of Candidate Points'])
disp(length(setCandidates))
disp(['setCandidates = '])
disp(setCandidates')


%%%%%%%%%%%%%%%%%5
%Creating the chart qPLAST-wECM,zECM
setPointsALL = cell(size(ECMdata_cluster)) ;

for icluster = 1:length(ECMdata_cluster)
    setPointsALL{icluster} =  ECMdata_cluster{icluster}.setPoints  ;
end
setPointsALL  =unique(cell2mat(setPointsALL' )) ;

ncluster = length(ECMdata_cluster) ;
npointsALL = length(setPointsALL) ;
wALL = zeros(npointsALL,ncluster) ;

for icluster = 1:ncluster
    setPloc = ECMdata_cluster{icluster}.setPoints  ;
    [dummy1,III,JJJ] = intersect(setPloc,setPointsALL,'stable') ;
    wALL(JJJ,icluster) =  ECMdata_cluster{icluster}.wRED  ;
end

ECMdata.setPoints = setPointsALL ;
ECMdata.wRED.Values = wALL ;
ECMdata.wRED.q = qPLAST ;
setElements_cand = large2smallINCLUDEREP(setPointsALL,DATA.MESH.ngaus) ;
ECMdata.setElements = setElements_cand ;

ECMdata.wRED.IndexDOFl_q = 2;



figure(4532)
hold on
title('WEIGHTS VERSUS LATENT VARIABLE (elastoplastic PROBLEM)')
xlabel('q')
ylabel('w')
for iddd  = 1:size(wALL,1)
    eeee = large2smallINCLUDEREP(setPointsALL(iddd),DATA.MESH.ngaus) ; 
    plot(qPLAST,wALL(iddd,:),'DisplayName',['Point =',num2str(setPointsALL(iddd)),' Elem = ',num2str(eeee)])
end
legend show



% 
% warning('BORRAR ESTO, ESTOY ENGAÑANDO AL PROGRAMA')
%     load('tmpWORKecm.mat','setPoints','wRED')
% ECMdata.setPoints = setPoints ;
% ECMdata.wRED.Values = repmat(wRED,1,length(qPLAST)) ;