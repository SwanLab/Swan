function [ECMdata_cluster,setCandidates] = ECM_local_clusteringTAILORED(nmodes,NEW_ORDER_MODES,wSTs,DATAoffline,BasisFint_permode,...
    DATA)

if nargin == 0
    load('tmp.mat')
end

setCandidates = [] ;
ECMdata_cluster = cell(1,nmodes) ;


for iclusterLOC = 1:nmodes
    icluster = NEW_ORDER_MODES(iclusterLOC) ;
    %     if icluster == 553
    %         disp('borrar esto....')
    %     end
    
    disp('*************************************************++')
    disp(['disp. mode = ',num2str(icluster)])
    disp('*************************************************++')
   
    [ECMdata_cluster{icluster},setCandidates ]=...
        ECM_clustersTAILORED(wSTs,DATAoffline,BasisFint_permode{icluster},DATA,setCandidates) ;
    disp(['*****************************'])
    
    
end