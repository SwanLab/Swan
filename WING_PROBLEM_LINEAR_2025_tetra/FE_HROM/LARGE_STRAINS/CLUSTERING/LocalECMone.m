function [ECMdata_cluster,setCandidates] = LocalECMone(NEW_ORDER_clusters,wFE,DATAoffline,A)

if nargin == 0
    load('tmp.mat')
end

nclusters = length(A) ; 
setCandidates = [] ;
ECMdata_cluster = cell(1,nclusters) ;
DATA = [] ; 


for iclusterLOC = 1:nclusters
    icluster = NEW_ORDER_clusters(iclusterLOC) ;
    %     if icluster == 553
    %         disp('borrar esto....')
    %     end
    
    disp('*************************************************++')
    disp(['CLUSTER = ',num2str(icluster)])
    disp('*************************************************++')
       [ECMdata_cluster{icluster},setCandidates ]=...
        LocalECMalg(wFE,DATAoffline,A{icluster},DATA,setCandidates) ;
    disp(['*****************************'])
    
    
end