function wDIAG = CompWeightDiag(wSTs,ndim)

wST = RepeatWEIGHT(wSTs,ndim); 
% Diagonal matrix with weights 
mTOT = length(wST) ; 
wDIAG = sparse(1:mTOT,1:mTOT,wST,mTOT,mTOT,mTOT) ;