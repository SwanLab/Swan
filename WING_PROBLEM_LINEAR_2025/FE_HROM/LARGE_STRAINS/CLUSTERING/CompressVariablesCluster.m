function  BasisU_cl_compress =  CompressVariablesCluster(BasisU_cl)

if nargin == 0
    load('tmp.mat')
end

% Compressing BasisU_r
TOL_BLOCK = zeros(length(BasisU_cl),1)' ;
DATAsvd=[];
[U,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(BasisU_cl,TOL_BLOCK,DATAsvd) ;

coeff = bsxfun(@times,V',S) ; 
 
