function [ECMdata_cluster,setCandidates,LIST_OF_CANDIDATES] = ...
    SAW_ECM(wFE,DATAoffline,BasisFint_cluster,DATA)

if nargin == 0
    load('tmp2.mat')
   % DATAoffline.NITERATIONS_NO_MORE_POINTS_negative_iterationECM = 10; 
  %  DATAoffline.TOL_ECM = 0.01; 
 end
nclusters = length(BasisFint_cluster) ; 
DATAoffline = DefaultField(DATAoffline,'SetCandidates_given_by_the_userECM',[]) ; 
setCandidates = DATAoffline.SetCandidates_given_by_the_userECM(:) ;
ECMdata_cluster = cell(1,nclusters) ;

AAA = tic ; 

  NEW_ORDER_clusters = nclusters:-1:1 ; 
LIST_OF_CANDIDATES = cell(1,nclusters) ; 
for iclusterLOC = 1:nclusters
    LIST_OF_CANDIDATES{iclusterLOC} = setCandidates ; 
    icluster =  NEW_ORDER_clusters(iclusterLOC) ;
    %     if icluster == 553
    %         disp('borrar esto....')
    %     end
    
    disp('*************************************************++')
    disp(['CLUSTER = ',num2str(icluster)])
    disp('*************************************************++')
     
    [ECMdata_cluster{icluster},setCandidates ]=...
        ECM_clustersNCns(wFE,DATAoffline,BasisFint_cluster{icluster},DATA,setCandidates) ;
    disp(['*****************************'])

     
end

AAA= toc(AAA) ; 
disp(['Time for the SAW-ECM: ',num2str(AAA),' s'])
