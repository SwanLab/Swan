function [  wSTs wST posgp] = ComputeW_RHS_cable(COOR,CN,TypeElement,nstrain,TypeIntegrand)
%%%%
% This subroutine   returns
%  1) wSTs = zeros(nelem*ngaus,1) ;  % Vector containinig the product of weights and Jacobians at all gauss points
%  2) wST = zeros(nelem*ngaus*nstrain,1) ; % Similar to the WeightST, but replicating   each entry nstrain times
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE),
% TypeElement: Type of finite element (quadrilateral,...),
% TypeIntegrand = 'K' or 'RHS'
%% Vectorized version
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 21-Apr-2016
%dbstop('10')
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
%TypeIntegrand = 'K';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
ngaus = length(weig) ;
% Initialization
%Belem = zeros(nstrain*nelem*ngaus,nnodeE*ndim) ;
wSTs = zeros(nelem*ngaus,1) ;  % Vector containinig the product of weights and Jacobians at all gauss points
wST = zeros(nelem*ngaus*nstrain,1) ; % Similar to the WeightST, but replicating the each entry nstrain times
% COORDINATE MATRIX arranged in a nelem*ndim x nnodeE matrix
XeALL= COORallELEM(ndim,nelem,nnodeE,CN,COOR) ;
idimREFERENCE = 1:ndim:size(XeALL,1);  % Only interpolation along this variable
XeALL = XeALL(idimREFERENCE,:) ;
% Let us define a matrix ROWSgauss such that   Belem(ROWSgauss(g,:),:) returns
% the B-matrices of the g-th points of all elements
indREF = 1:nstrain*ngaus*nelem ;
ROWSgauss = reshape(indREF,nstrain,nelem*ngaus) ;
weigREP = repmat(weig',nelem,1)  ; % nelem x 1 tiling copies of weig
for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ;
    % Jacobian Matrix for the g-th G. point of all elements %
    JeALL = XeALL*BeXi' ;
    %%%%%%%%%
    % JAcobian
    detJeALL= JeALL ; 
    % Matrix of derivatives with respect to physical coordinates
    %   inv_JeTall = inverseTRANSvectorize(JeALL,ndim,detJeALL) ;
    %   BeTILDEall = inv_JeTall*BeXi ;
    % Matrix of symmetric gradient
    % dbstop('18')
    %  BeALL = QtransfBvect(BeTILDEall,ndim) ; % Transformation from B-matrix for scalar fields to B-matrix for vector fields
    % Assigning BeALL ( g-th Gauss point) to the global matrix Belem
    ROWSglo =  ROWSgauss(:,g:ngaus:ngaus*nelem);
    ROWSglo = ROWSglo(:) ;
    % Belem(ROWSglo,:) = BeALL ;
    % Weight vectors
    % --------------
    [wST,wSTs] = DetermineWeightsST(detJeALL,weigREP,ngaus,g,nelem,wST,wSTs,nstrain,ROWSglo);
    
end
