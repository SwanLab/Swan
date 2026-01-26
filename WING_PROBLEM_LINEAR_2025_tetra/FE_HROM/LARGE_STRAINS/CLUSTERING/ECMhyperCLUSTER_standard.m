function [ECMdata_cluster,TableInfoModesCluster] = ECMhyperCLUSTER_standard(BasisU_cluster,BasisPone_cluster,BasisStwo_cluster,...
    OPERFE,DISP_CONDITIONS,DATA_GENGAUSS,DATA,DATAoffline)

nclusters = length(BasisU_cluster) ; 

ECMdata_cluster = cell(1,nclusters) ; 
%DATAoffline.TryUseSameECMpointsForAllClusters  == 1

setAllPoints = cell(1,nclusters) ; 

for icluster = 1:nclusters
    disp(['CLUSTER = ',num2str(icluster)])
    
    BasisU = BasisU_cluster{icluster} ;
    BasisPone = BasisPone_cluster{icluster} ;
    BasisStwo = BasisStwo_cluster{icluster} ;
    
   ECMdata_cluster{icluster} = ECM_clusters(OPERFE,DISP_CONDITIONS,BasisU,DATAoffline,DATA_GENGAUSS,BasisPone,DATA,BasisStwo) ; 
    
        disp(['*****************************'])
        
    setAllPoints{icluster} = ECMdata_cluster{icluster}.setPoints' ;     
end

setAllPoints_unique = unique(cell2mat(setAllPoints)) ; 

disp('********************************************************************')
disp(['Total number of points   =',num2str(length(setAllPoints_unique)),' (out of ',num2str(length(cell2mat(setAllPoints))),')']) ;
disp('********************************************************************')

