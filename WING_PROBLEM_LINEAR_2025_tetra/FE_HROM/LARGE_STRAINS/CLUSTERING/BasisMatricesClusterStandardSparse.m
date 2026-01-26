function [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline] = BasisMatricesClusterStandardSparse(DATAoffline,SNAPdisp,DISP_CONDITIONS)

if nargin == 0
    load('tmp1.mat')
end
DISTANCE_CENT = [] ;  DISTANCES = [] ; 
CENT = [] ; 
    INDEX_SELECT  = 1:size(SNAPdisp,2) ; 
    
           % Remove zero points 
           NNN = sum(SNAPdisp.^2,1) ; 
          ddddd =  find(NNN<1e-16) ; 
          INDEX_SELECT(ddddd) = [] ;  
      %    SNAPdisp(:,ddddd) = [] ; 
          NNNdiff = diff(ddddd) ; % Provisional, modify this 
          DATAoffline.InitialclusterINdex = [2:NNNdiff(1):size(SNAPdisp,2)] ; 

[CCC,INDcluster] = SSC_jaho(SNAPdisp(:,INDEX_SELECT),DATAoffline.NCLUSTERS_BASIS_DISP,DATAoffline) ;

  IND =  ones(size(SNAPdisp,2),1) ;   
          IND(INDEX_SELECT) = INDcluster ; 



DATAoffline.InitialclusterINdex = IND(DATAoffline.InitialclusterINdex )   ; %%%%s

% LOOP  OVER CLUSTERS TO COMPUTE the  BASIS MATRICES
% --------------------------------------------

nclusters = DATAoffline.NCLUSTERS_BASIS_DISP; 
BasisU_cluster = cell(1,nclusters) ;
errorCLUSTER = zeros(1,nclusters) ;
% DISTANCES = zeros(nclusters) ;

WITH_ALL_DOFS = DATAoffline.CLUSTER.WITH_ALL_DOFS  ;
%
if  WITH_ALL_DOFS == 1
    DOFSincluded = 1:size(SNAPdisp,1) ;
else
    DOFSincluded = DISP_CONDITIONS.DOFl ;
    
end

for i=1:nclusters
    
    STEPSCLUSTER = find(IND == i) ;
    % [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA) ;
    e0 = DATAoffline.errorDISP_cluster ;
    mu = 0 ; R = 0 ; DATALOC.RELATIVE_SVD = 1;
    DATALOC.HIDE_OUTPUT = 1 ;
    [BasisU_cluster{i},SS,VV,errorCLUSTER(i)] =   RSVDT(SNAPdisp(DOFSincluded,STEPSCLUSTER),e0,mu,R,DATALOC) ;
    
    disp(['SVD cluster =',num2str(i),'; Nmodes =',num2str(length(SS))]) ;
%     
%     for j = 1:nclusters
%         DISTANCES(i,j) = norm(CENT(i,:)-CENT(j,:));
%     end
    
    
end

DATAoffline =  DefaultField(DATAoffline,'ExploreIsoMap',0) ;


if DATAoffline.ExploreIsoMap == 1
    
    ExploreIsoMapForClustering(SNAPdisp,DATAoffline) ;
end