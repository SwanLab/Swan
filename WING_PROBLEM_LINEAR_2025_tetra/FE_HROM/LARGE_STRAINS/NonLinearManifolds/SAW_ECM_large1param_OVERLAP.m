function ECMdata = SAW_ECM_large1param_OVERLAP(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT)


 
SNAPfint = cell2mat(SNAPfint) ; 
[qLATENT,ii] = sort(qLATENT) ; 
[aab,bbb] = find(abs(qLATENT)<1e-10) ; 
if ~isempty(bbb)
    qLATENT(bbb) = [] ; 
    ii(bbb) = [] ; 
end


SNAPfint = SNAPfint(:,ii) ; 
DATAoffline = DefaultField(DATAoffline,'ECM_Number_Snapshots_Overlapping',1) ;
DATAoffline = DefaultField(DATAoffline,'IntegrateExactlyVolume_SAW_ECM',1) ;

SNAPfint_w  = cell(1,size(SNAPfint,2)) ; 
 
sqW = sqrt(OPERFE.wSTs) ; 
sqWfe = diag(sparse(sqW)) ; 
NumberModesCluster = zeros(size(SNAPfint_w)) ; 
 INCLUDE_sum_vol  =DATAoffline.IntegrateExactlyVolume_SAW_ECM; 
for icluster = 1:size(SNAPfint,2)
 
    disp(['icluster = ',num2str(icluster)])
%     if icluster == 299
%         disp('')
%     end
    icluster_back = max(1,icluster-DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluster_forw = min(size(SNAPfint,2),icluster+DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluster_select = [icluster_back:icluster_forw] ;
     
        SNAPlocI =  SNAPfint(:,icluster_select) ;
        % INCLUDE ELASTIC SNAPSHOTS/integ. exact volum
        SNAPlocI = sqWfe*SNAPlocI ;  
        if INCLUDE_sum_vol ==0
                      
            TOL = [DATAoffline.errorFINT] ;
            DDD.ISRELATIVE = 1;
            [UU,SS,VV] = SVDT(SNAPlocI,TOL,DDD) ;
        else
              TOL = [0,DATAoffline.errorFINT] ;
              DDD.HIDE_OUTPUT = 1; 
            [UU,SS,VV] = SRSVD({sqW,SNAPlocI},TOL,DDD) ;
        end
    
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
%Creating the chart qLATENT-wECM,zECM
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

figure(4532)
hold on
title('WEIGHTS VERSUS LATENT VARIABLE (HYPERELASTIC PROBLEM)')
xlabel('q')
ylabel('w')
for iddd  = 1:size(wALL,1)
    eeee = large2smallINCLUDEREP(setPointsALL(iddd),DATA.MESH.ngaus) ; 
    plot(qLATENT,wALL(iddd,:),'DisplayName',['Point =',num2str(setPointsALL(iddd)),' Elem = ',num2str(eeee)])
end
legend show

ECMdata.setPoints = setPointsALL ;
ECMdata.wRED.Values = wALL ;
ECMdata.wRED.q = qLATENT ;
setElements_cand = large2smallINCLUDEREP(setPointsALL,DATA.MESH.ngaus) ;
ECMdata.setElements = setElements_cand ;

ECMdata.wRED.IndexDOFl_q = 2;

% 
% warning('BORRAR ESTO, ESTOY ENGAÑANDO AL PROGRAMA')
%     load('tmpWORKecm.mat','setPoints','wRED')
% ECMdata.setPoints = setPoints ;
% ECMdata.wRED.Values = repmat(wRED,1,length(qLATENT)) ;