function [InvCNmatrix ElemNode maxELEM]= CondenseSparseMatrix(InvCN,nelem,nnode)


% Number of element per node
ElemNode = full(sum(InvCN,2)) ;
% Maximum number of elements per nod
maxELEM = full(max(ElemNode)) ;
% Nonzero elements
[NonZeroElem] = find(InvCN') ;
[NonZeroElemColum ROW ] = ind2sub([nelem,nnode],NonZeroElem);
% Sought matrix
InvCNmatrix= zeros(nnode,maxELEM,1);
ElemNodeIni = [cumsum(ElemNode)] ;
 ElemNodeIni = [1+ElemNodeIni(1:end-1),];
 ElemNodeIni =[1;ElemNodeIni];
for ielem = 1:maxELEM
    % ielem associated to each nnode
    % First we select the elements that have equal or more than ielem for
    % node
    IndElemLoc = find(ElemNode>=ielem) ;
    % In the set of nonzero entries of InvCN, we now determine the 
    % local numbering through  ElemNodeIni
    IndElemGLO = ElemNodeIni(IndElemLoc)+ielem-1 ; 
    % Thus
    InvCNmatrix(IndElemLoc,ielem) = NonZeroElemColum(IndElemGLO) ;
    
end