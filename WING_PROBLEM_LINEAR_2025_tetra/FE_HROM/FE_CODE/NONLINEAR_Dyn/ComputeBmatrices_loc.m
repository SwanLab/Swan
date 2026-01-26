function [ wSTs  XeALL  Bst wST  posgp]= ...
    ComputeBmatrices_loc(COOR,CN,TypeElement,DATA)
%%%%
%  Copy of function devoted to generate STiffness matrix in the elastic
%  regime . JAHO, 20-sept-2018
% ----------------
% This subroutine   returns:
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
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015, JAHO, 20-sept-2018
%dbstop('11')
if nargin == 0
    load('tmp.mat')
end


%dbstop('36')
nstrain = DATA.nstrain ;
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
%------------------------------------------------
% Matrix consisting of all element B-matrices (wST -->Weights)
time1 = tic ;
disp('Computing B-matrices for all elements...')
time2 = tic ;
%dbstop('44')
[ Belem wSTs wST XeALL posgp] = ComputeBelemALL(COOR,CN,TypeElement,nstrain,DATA) ;
% -------------------------------------------------------------------------
time2 = toc(time2) ;
disp(['Done (in =',num2str(time2),' s)']);
% ---------------------------------------------------------------------------
ngaus = length(wSTs)/nelem ; 

% Global matrix such that STRAIN_aLL_gauss_points = Bst*d 
Bst = AssemblyBGlobal(Belem,nstrain,nelem,nnodeE,ndim,ngaus,CN,nnode) ;
% 
% 
% [K,Bst,Cglo] = AssemblyMethodBCB(wSTs,Belem,nstrain,nelem,nnodeE,ndim,CN,nnode,celasglo,wST) ;
% %
% % else
%     % Standardard, vectorized
%     %  BstW = [] ; Bst = [] ;
%     [K,Bst,Cglo]= AssemblyMethodStandard(wSTs,Belem,nstrain,nelem,nnodeE,ndim,CN,nnode,celasglo,wST);
% end
% %timeK = toc(time1);
% %dbstop('55')
% %disp(['Elapsed time = ',num2str(timeK),' s, for BCBformulation = ',num2str(DATA.BCBformulation)])
%
%
%
%
%
%
%
