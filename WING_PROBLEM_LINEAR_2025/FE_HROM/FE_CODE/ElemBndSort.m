function  setBelem = ElemBndSort(CONNECTb,CNb)
% CONNECTb=  Connectivity matrix for   boundary elements
% CNb =  Connectivity matrix for  a subset of  boundary elements 
% setBelem --> Numbering of CNb within CONNECTb
% Modified (or at least revisited, 13 Dec. 2016,  )
% There is a vectoriced version: ElemBnd. ... NOT TRUE (23-Aug-2019)
%
if nargin == 0
    load('tmp1.mat')
end
setBelem = zeros(size(CNb,1),1) ; 
for i=1:size(CNb,1)
    ComparD =  bsxfun(@minus,CONNECTb,CNb(i,:)) ;      
    a = find(sum(abs(ComparD),2)==0) ;  % I forgot the "abs" ... : ), 13-Dec-2016
    setBelem(i) = a(1) ; 
end