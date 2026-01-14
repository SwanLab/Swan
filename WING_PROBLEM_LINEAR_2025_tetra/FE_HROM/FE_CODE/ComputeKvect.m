function [K CN wSTs  XeALL Cglo Bst wST MaterialType IndicesRenumberingElements]= ...
    ComputeKvect(COOR,CN,TypeElement, celasglo,DATA,MaterialType)
%%%%
% This subroutine   returns:
% 1) the global stiffness matrix K (ndim*nnode x ndim*nnode)
% 2)  wSTs = zeros(nelem*ngaus,1) ;  % Vector containinig the product of weights and Jacobians at all gauss points
%  4)XeALL:  COORDINATES of all elements, arrangedin a nelem*ndim x nnodeE matrix
%  5) CN: Connectivity matrix (renumbered so as to minimize the bandwidth of Bst)
%  6) Cglo: Blockwise diagonal matrix containing the elasticity matrices of
%  all gauss points
%  7) Bst: Sparse matrix such that strainST = Bst*d, where strainST = [strain(e=1,g=1); strain(e=1,g=2)...
%                                                                   strain(e=1,g=1); strain(e=2,g=2) ...]
% --------------------------------------------------------------------------------------
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE),
% TypeElement: Type of finite element (quadrilateral,...),  celasglo (nstrain x nstrain x nelem)
% Array of elasticity matrices
% % 6. Miscellaneous input DATA  -->
% DATA.BCBformulation = 1;  % Method for assembly K . If == 1, we make K =
% BstW^T*Cglo*BstW. Otherwise, K is assembled in the standard fashion (yet in a vectorized fashon)

%% Vectorized version
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015
%dbstop('11')
if nargin == 0
    load('tmp1.mat')
end



if  DATA.RENUMBERED_OUTSIDE == 0
    [~,IndicesRenumberingElements]  = sort(CN(:,1)) ;
    CN = CN(IndicesRenumberingElements,:) ;
    celasglo = celasglo(:,:,IndicesRenumberingElements) ;
    if ~isempty(MaterialType)
        MaterialType = MaterialType(IndicesRenumberingElements) ;
    end
else
    IndicesRenumberingElements = [] ;
end



%error('Change thissssssssssssss .... DOmAINVAR')
%%% Dimensions of the problem
%dbstop('36')
if ~isstruct(celasglo)
    nstrain = size(celasglo,1) ;
else
    nstrain = size(celasglo.Celas,2) ;
    
end

nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
%------------------------------------------------
% Matrix consisting of all element B-matrices (wST -->Weights)
time1 = tic ;
disp('Computing B-matrices for all elements...')
time2 = tic ;
%dbstop('44')
[ Belem wSTs wST XeALL] = ComputeBelemALL(COOR,CN,TypeElement,nstrain,DATA) ;
% -------------------------------------------------------------------------
time2 = toc(time2) ;
disp(['Done (in =',num2str(time2),' s)']);

DATA = DefaultField(DATA,'BCBformulation',1) ; 

if DATA.BCBformulation == 1
    % Method in which the stiffness matris is computed as: K =% BstW^T*Cglo*BstW.
    % ---------------------------------------------------------------------------
    
    [K,Bst,Cglo] = AssemblyMethodBCB(wSTs,Belem,nstrain,nelem,nnodeE,ndim,CN,nnode,celasglo,wST) ;
    
else
    % Standardard, vectorized
    %  BstW = [] ; Bst = [] ;
    [K,Bst,Cglo]= AssemblyMethodStandard(wSTs,Belem,nstrain,nelem,nnodeE,ndim,CN,nnode,celasglo,wST);
end
%timeK = toc(time1);
%dbstop('55')
%disp(['Elapsed time = ',num2str(timeK),' s, for BCBformulation = ',num2str(DATA.BCBformulation)])







