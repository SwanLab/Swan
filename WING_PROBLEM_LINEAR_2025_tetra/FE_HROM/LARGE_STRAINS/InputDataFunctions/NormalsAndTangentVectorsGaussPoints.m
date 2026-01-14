function [UnitNormals,tau] = NormalsAndTangentVectorsGaussPoints(MESH,CNb)
%--------------------------------------------------------------------------
%  NormalsAndTangentVectorsGaussPoints
%
%  Computes the **tangent** and **normal unit vectors** at each Gauss point
%  of the boundary elements defined by the connectivity matrix `CNb`.
%  These vectors are required for the evaluation of boundary integrals,
%  including:
%    - Neumann boundary conditions (traction vectors)
%    - Fluid–structure coupling (normal projections)
%    - Follower loads and hydrostatic pressure contributions
%
%  INPUTS:
%    - MESH : structure containing:
%        > COOR          : nodal coordinates
%        > TypeElementB  : type of boundary element
%    - CNb  : connectivity matrix for boundary elements
%
%  OUTPUTS:
%    - UnitNormals : unit outward normal vector at each Gauss point
%    - tau         : tangent vectors at each Gauss point
%                    (1 in 2D, 2 in 3D as a cell array {τ₁, τ₂})
%
%  STEPS:
%    1. Evaluate shape functions and derivatives at boundary Gauss points.
%    2. Build `Lbool`, the Boolean operator mapping nodal to Gauss point coordinates.
%    3. Build directional transformation matrices `BstB` using `dershapef`.
%    4. Compute tangent vectors at Gauss points using:
%           τᵢ = BstB{i} * Lbool * x
%    5. Compute unit normals as:
%         - 2D:  n = [τ₂ ; -τ₁]
%         - 3D:  n = normalize(cross(τ₁, τ₂))
%
%  REMARKS:
%    - Assumes **all boundary elements belong to the same family** (same interpolation).
%    - Gauss point indexing is handled in block format using `ConvertBlockDiag_general`.
%    - Vector `X = COOR'` is used to construct the global geometry vector.
%
%  APPLICATIONS:
%    - Traction vector projection onto normals: t(x_gp) · n
%    - Hydrostatic and pressure follower loads
%    - Weak enforcement of boundary conditions in variational forms
%
%  AUTHOR:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 31-Aug-2021
%    Comments by ChatGPT4, 31-May-2025
%  SEE ALSO:
%    - ComputeElementShapeFun
%    - ConvertBlockDiag_general
%    - Lbool_vectorized
%
%--------------------------------------------------------------------------




% Computation of tangent and normal vector (unitary) at each Gauss point
% of the surface defined by the table of connectivities CNb
% JAHO, 31-Aug-2021,


nelemB = size(CNb,1) ; % Number elementes boundary
nnodeE = size(CNb,2) ; % Number of nodes per element (boundary)

TypeIntegrand = 'RHS';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(MESH.TypeElementB,nnodeE,TypeIntegrand) ;

ngaus = size(posgp,2) ;   % Number of Gauss points.
nnode = size(MESH.COOR,1) ; % Number of nodes, for all the discretization
ndim = size(MESH.COOR,2);
Lbool = Lbool_vectorized(CNb,nnode,ngaus,ndim)  ;
%Lbool = Lbool_nonvectorized(CNb,nnode,ngaus,ndim)  ;

DATA.MESH.HYDRO.nelemB = nelemB ; % Number of elements wet surface
DATA.MESH.HYDRO.ngausT = nelemB*ngaus ; % Number of elements wet surface
DATA.MESH.HYDRO.ngaus = ngaus ;
DATA.MESH.HYDRO.nnodeE = nnodeE ;
% ------------------------------------------------------
% 2. MAtrix of shape functions for ALL Gauss points (concatenation, and multiplied by the weights)
% -------------------------------------------------
%wSTft =repmat(weig',nelemB,1) ;      ; % Vector of weights for all Gauss points


% Nescl_w = bsxfun(@times,shapef,weig') ;  % Scalar functions multiplied by the weights
%
% % Each entry multiplied by the identiy  N1*IDENT, N2*IDENT2 ....
% Nb_xi_w = zeros(ndim*size(Nescl_w,1),ndim*size(Nescl_w,2)) ;
% for igaus = 1:size(Nescl_w,1)
%     INI = (igaus-1)*ndim+1; FIN = igaus*ndim ;
%     Nb_xi_w(INI:FIN,:) = StransfN(Nescl_w(igaus,:),ndim) ;
% end

% Since we are assuming that all boundary elements are of the same family,
% this matrix consists in NGAUS_all tiling copies of shapef
% Strictly speaking, this would   not be necessary  in this
% case, as we are assuming that all elements pertains to the same family (we do it for the sake
% of tenerality)
%NbST_w_elem = repmat(Nb_xi_w, nelemB,1) ;



% ---------------------------------
% 3. Convert NbST_w_elem  into a diagonal matrix
%
% Now we have NbST_w_elem = [N_1; N_2 ... N_M]  (M gauss points), and we
% wish to turn it into a block diagonal matrix
% We have to adapt /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/ConvertBlockDiag.m% to cope
%   *** /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/ConvertBlockDiag_general.m
% [NbST_diag_w,irowsNDIM,icolsNDIM ]= ConvertBlockDiag_general(NbST_w_elem,ndim) ;
%
% NbST_w =  NbST_diag_w*Lbool ;


% 4.  Matrix BstB_allgauss, such that   tau_i = diag(BstB_allgauss{i})*Lbool*COOR'(:) gives the
% tangent vectors in the i-direction
% \BstB{i}{e}
% --------------------------------------------------------------------------
% For a given element
% -------------------
BstB_e = cell(ndim-1,1) ;
BstB_e{1} = zeros(ndim*size(shapef,1),ndim*size(shapef,2)) ;
for igaus = 1:size(shapef,1)
    INI = (igaus-1)*ndim+1; FIN = igaus*ndim ;
    for  idim  = 1:ndim-1
        BstB_e{idim}(INI:FIN,:) = StransfN(squeeze(dershapef(idim,:,igaus)),ndim) ;
    end
end
% Since all elements are of the same family
BstB_allgauss = cell(ndim-1,1) ;
for  idim  = 1:ndim-1
    BstB_allgauss{idim}  = repmat(BstB_e{idim}, nelemB,1) ;
end

% -------------------------------------------------------------------------
% 5. Initial tangent vectors
% -------------------------
tTANGiniST = cell(1,ndim-1) ;
X = MESH.COOR' ;
for idim = 1:ndim-1
    tTANGiniST{idim} = ConvertBlockDiag_general(BstB_allgauss{idim},ndim)*Lbool*X(:) ;
end
% See assessment in FLOAT_IMPLE.pdf, kw:89


% 6. Initial normals
% ------------------
tau1 = reshape(tTANGiniST{1},ndim,[]) ;
tau1_n  = sqrt(sum(tau1.^2,1)) ;
tau1 = bsxfun(@times,tau1',1./tau1_n')';

if ndim == 3
    tau2 = reshape(tTANGiniST{2},ndim,[]) ;
    
    % tau2_n  = sqrt(sum(tau2.^2,1)) ;
    % tau2 = bsxfun(@times,tau2',1./tau2_n')';
    
    
    mINI =  cross(tau1,tau2) ;
    
    nmINI  = sqrt(sum(mINI.^2,1)) ;
    UnitNormals = bsxfun(@times,mINI',1./nmINI')';
    
    tau2 = cross(UnitNormals,tau1) ;
    
    tau = {tau1,tau2};
    
else
    tau = {tau1} ; 
    
    UnitNormals  = [tau1(2,:)
                    -tau1(1,:)] ; 
    
end
