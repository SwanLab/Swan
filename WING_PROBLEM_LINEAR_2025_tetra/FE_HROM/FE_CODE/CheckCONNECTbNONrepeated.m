function [CONNECTbNonRepeated,NODES_repeated]= CheckCONNECTbNONrepeated(CONNECTb)
% Given a connectivity matrix CONNECTb, CheckCONNECTbNONrepeated
% returns the non-repeated elements of the matrix 
% Created JAN-7th-2019 to amend an issue of GID. 
% It turns out that when creating the Boundary Elements, Gid sometimes
% provides repeated elements 

if nargin == 0
    load('tmp.mat')
  %  CONNECTb = CONNECTbNonRepeated ; 
end

[CONNECTbNEW, IXROWS ]= sort(CONNECTb,2) ; 

[AA,IDX] = sort(CONNECTbNEW(:,1)) ; 

CONNECTbNEW = CONNECTbNEW(IDX,:) ; 

LLL = sum(abs(diff(CONNECTbNEW,1)),2) ;
LLL = find(LLL==0) ; 
% Elements to remove 
NODES_repeated = [] ; 
if ~isempty(LLL)
    ElemsToRemove  = IDX(LLL) ;
    
    NODES_repeated= unique(CONNECTbNEW(LLL,:)); 
    LIST = 1:size(CONNECTb,1) ; 
    LIST(ElemsToRemove) = [] ; 
    
    CONNECTbNonRepeated = CONNECTb(LIST,:) ; 
else
    CONNECTbNonRepeated = CONNECTb ;
end