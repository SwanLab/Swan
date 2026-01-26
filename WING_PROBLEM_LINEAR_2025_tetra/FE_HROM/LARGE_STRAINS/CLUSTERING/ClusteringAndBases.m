function [BasisU_cluster,INDICES_CLUSTERS,IND,CentroidIndex] = ClusteringAndBases(A,DATAoffline)
 DATAoffline =  DefaultField(DATAoffline,'DistanceClustering','sqeuclidean') ;


 DATAoffline = DefaultField(DATAoffline,'USE_INITIAL_STARTING_POINT_CLUSTERING',1) ; 
 
 if DATAoffline.USE_INITIAL_STARTING_POINT_CLUSTERING == 1
     INDICES_snapshots = DEIMfun_SELECT(A,DATAoffline) ;
     [IND,CENT,~,DISTANCE_CENT] = kmeans(A',DATAoffline.NCLUSTERS_BASIS_DISP,...
         'start',A(:,INDICES_snapshots)','Distance',DATAoffline.DistanceClustering) ;
 else
     SomeParticularNumericSeed = 1; 
     rng(SomeParticularNumericSeed)
     [IND,CENT,~,DISTANCE_CENT] = kmeans(A',DATAoffline.NCLUSTERS_BASIS_DISP,'Distance',DATAoffline.DistanceClustering) ;
 end
 
 [~,CentroidIndex ]= min(DISTANCE_CENT,[],1) ; 
 


BasisU_cluster = cell(1,DATAoffline.NCLUSTERS_BASIS_DISP) ; 
INDICES_CLUSTERS = cell(1,DATAoffline.NCLUSTERS_BASIS_DISP) ;  
for i=1:size(CENT,1)    
    STEPSCLUSTER = find(IND == i) ;  INDICES_CLUSTERS{i} =  STEPSCLUSTER ; 
    % [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA) ;
    e0 = DATAoffline.TOL_SVD_MATRIX;
    mu = 0 ; R = 0 ; DATALOC.RELATIVE_SVD = 1;
    DATALOC.HIDE_OUTPUT = 1 ;
    [BasisU_cluster{i},SS,VV,errorCLUSTER] =   RSVDT(A(:,STEPSCLUSTER),e0,mu,R,DATALOC) ;
    disp(['SVD cluster =',num2str(i),'; Nmodes =',num2str(length(SS))]) ;   
end