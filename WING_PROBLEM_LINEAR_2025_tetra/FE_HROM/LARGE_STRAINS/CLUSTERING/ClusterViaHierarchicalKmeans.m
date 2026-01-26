function [BasisU_cluster,INDICES_CLUSTERScell,CentroidIndex,CentroidCoordinates]  = ClusterViaHierarchicalKmeans(A,DATAoffline,DISP_CONDITIONS)
% Hierarchical k-means 
% INPUTS 
% ********
% 1) A: Snapsthot matrix
% 2) DATAoffline.nminBASIS      Maximum number of displacement modes per
% cluster 
% 3) DATAoffline.NCLUSTERS_BASIS_DISP   :   Initial number of clusters  (Default =1)
% 4) DATAoffline.TOL_SVD_MATRIX      % Relative tolerance SVD  (Default = 0)
% ******************
% OUTPUTS 
%***********+
% BasisU_cluster: Basis matrices for each clusterb 
% INDICES_CLUSTERScell : Indexes (columns) for each cluster
% CentroidIndex : Indexes centroids 

%-----------------------------------
% Joaquin A. Hernandez, 9-Feb-2021 

% 
 
  


DATAoffline = DefaultField(DATAoffline,'NCLUSTERS_BASIS_DISP',1) ; 
DATAoffline = DefaultField(DATAoffline,'errorDISP_cluster',0) ; 



[BasisU_cluster,INDICES_CLUSTERScell,IND,CentroidIndex,CentroidCoordinates] = ClusteringAndBases(A,DATAoffline,DISP_CONDITIONS) ;
EXIT = 0  ; 
kiter = 1; 
while  EXIT == 0 
    disp('************************************************')
    disp(['REFINEMENT ITERATION = ',num2str(kiter)] )
     [BasisU_cluster,INDICES_CLUSTERScell,CentroidIndex,CentroidCoordinates] = ClusterRefinement(A,BasisU_cluster,INDICES_CLUSTERScell,DATAoffline,CentroidIndex,CentroidCoordinates,DISP_CONDITIONS) ; 
     [AA,nMODESall] = cellfun(@size,BasisU_cluster) ;     
     if all(nMODESall<= DATAoffline.nminBASIS)
         EXIT = 1; 
     end         
      disp('************************************************')
end