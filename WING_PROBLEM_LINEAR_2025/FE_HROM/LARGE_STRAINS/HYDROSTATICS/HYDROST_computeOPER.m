function   [HYDRO,DATA,MESH] = HYDROST_computeOPER(MESH,DATA)
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Build the OFFLINE operators required to compute hydrostatic forces on a
%   floating body. The routine assembles boundary Gauss-point data, tangent
%   vectors, normals, shape-function operators, boolean scatterers, and
%   weighted quadrature terms for subsequent ONLINE evaluations of pressure
%   and follower loads.
%
% INPUTS:
%   MESH : structure with discretization data (3D domain assumed here).
%          • MESH.COOR        → (nnode × 3) nodal coordinates.
%          • MESH.CNb         → (nelemB × nnodeEb) boundary connectivity.
%          • MESH.TypeElementB→ boundary element family (line/tri/quad on the surface).
%          • MESH.Indexes_faces_bnd_element (optional) to pick wet faces.
%   DATA : structure with problem data.
%          • DATA.FOLLOWER_LOADS.HYDROSTATIC.WET_SURFACES (optional list of labels).
%          • DATA.vGRAVITY     → gravity direction (expects X2-direction nonzero).
%
% OUTPUTS:
%   HYDRO : structure containing OFFLINE hydrostatic operators:
%           • HYDRO.tTANGiniST      → {1:2} cell; initial (non-unit) tangent vectors τ₁, τ₂
%                                     at each boundary Gauss point (includes Jacobian).
%           • HYDRO.mNORMiniST      → stacked (non-unit) normals m = τ₁ × τ₂ at Gauss points.
%           • HYDRO.BstB            → {1:2} cell; derivative shape matrices along surface
%                                     parametric directions (per Gauss point, all elements).
%           • HYDRO.yPRESSiniST     → initial X2-coordinate at Gauss points (for hydrostatic head).
%           • HYDRO.Lbool           → boolean scatterer from nodal dofs to boundary Gauss points.
%           • HYDRO.NbSTgrav        → operator mapping nodal displacements to X2-increments
%                                     at Gauss points (gravity direction component of N).
%           • HYDRO.NbST_w          → weighted shape-function operator for pressure resultants
%                                     (so that F_pr = −NbST_wᵀ · t_pr).
%           • HYDRO.NbST            → unweighted boundary shape-function operator (vector field form).
%           • HYDRO.wST             → stacked Gauss weights (tiled per element).
%           • HYDRO.JacobianWeights → ||m|| = ||τ₁ × τ₂|| (surface Jacobian per Gauss point).
%           • HYDRO.irows1, icols1  → indices from block-diagonalization (scalar case).
%           • HYDRO.irowsNDIM, icolsNDIM → indices from block-diagonalization (vector case).
%   DATA  : updated defaults (e.g., WET_DAMPING_COEFFICIENT if absent).
%   MESH  : echoes HYDRO boundary connectivity in MESH.HYDRO.CNb.
%
% METHOD / PIPELINE (vectorized):
%   0) Select “wet” boundary faces: if DATAHYDRO.WET_SURFACES is empty, use all;
%      otherwise gather and stack the chosen face connectivities into CNb.
%   1) Quadrature & shapes on boundary elements (TypeIntegrand='RHS'):
%        [weig, posgp, shapef, dershapef] = ComputeElementShapeFun(...)
%   2) Build Gauss-point boolean scatterer:
%        Lbool = Lbool_vectorized(CNb, nnode, ngaus, ndim=3)
%   3) Weighted shape functions:
%        Nescl_w = shapef ⊗ I weighted by w_g  → NbST_w_elem (vector field form)
%        ConvertBlockDiag_general → NbST_diag_w; then NbST_w = NbST_diag_w * Lbool
%   4) Surface directional derivatives (per parametric dir i=1,2):
%        BstB_e{i} from dershapef(i,:); replicate for all elements → BstB_allgauss{i}
%   5) Initial tangents & normals:
%        tTANGiniST{i} = ConvertBlockDiag_general(BstB_allgauss{i},...) * Lbool * vec(COOR)
%        m = τ₁ × τ₂; store HYDRO.mNORMiniST and its norm as JacobianWeights
%   6) Gravity-direction projector (expects DATA.vGRAVITY(2) ≠ 0):
%        NbSTgrav extracts the X2 component of N at Gauss points; then
%        yPRESSiniST via block-diagonalization of NbSTgrav and application to COOR.
%   7) Also assemble the unweighted boundary operator NbST (vector field form) for clarity.
%
% ASSUMPTIONS / SANITY:
%   - 3D domain with 2D boundary elements (ndim=3, ndimB=2).
%   - All boundary elements share the same interpolation family (tiling shapes is valid).
%   - Gravity must have a nonzero X2 component; otherwise an error is raised.
%   - Normals and tangents are “non-unit” here; magnitudes include the surface Jacobian.
%
% FILE REFERENCES (internal notes):
%   - DynamicFloatingBodies.tex
%   - MLEARN_develop.pdf
%   - FLOAT_IMPLE.pdf
%
% AUTHOR / HISTORY:
%   Initial implementation: 1-Jul-2021
%   Comments clarification: 7-Nov-2025
%   JAHO — Joaquín A. Hernández — jhortega@cimne.upc.edu
% -------------------------------------------------------------------------




% Function for computing the OFFLINE operators required for
% determining hydrostatic forces
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
% PDF:
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/MLEARN_develop.pdf
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/FLOAT_IMPLE.pdf
%
%  
% HYDRO.tTANGiniST = tTANGiniST ;   % Tangent vectors at each Gauss point (not unitary)
% HYDRO.mNORMiniST = mNORMiniST ;   % Normal vector at each GAuss point (not unitary, it includes the Jacobian of the Gauss point )

% HYDRO.BstB = BstB_allgauss;  % Derivate shape function at all Boundary Gauss points   

% HYDRO.yPRESSiniST = yPRESSiniST;  % Initial coordinate X2 
% HYDRO.Lbool = Lbool;  % Boolean operator 
% HYDRO.NbSTgrav = NbSTgrav;  %  MAtrix such that NbSTgrav*Lbool*d gives the increment of x2
% HYDRO.NbST_w = NbST_w;  %  MAtrix such that Fpr = -NbST_w^t*t_pr   
%
%
% JAHO, 1-July-2021
% -------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
DATAHYDRO = DATA.FOLLOWER_LOADS.HYDROSTATIC ;

% Connectivity matrix (boundary)
%CNb_all = MESH.CNb ;


% 0. Determine which are the (potentially) "wet" surfaces, and construct
% the corresponding matrix of connectivities
% ---------------------------------------------------------

DATAHYDRO = DefaultField(DATAHYDRO,'WET_SURFACES',[]) ;

if isempty(DATAHYDRO.WET_SURFACES)
    CNb = MESH.CNb ;
    MaterialLibraryBND = ones(size(CNb,1)) ;
else
    IndexesWetSurfaces = DATAHYDRO.WET_SURFACES ;
    IndexElemSelect = MESH.Indexes_faces_bnd_element(IndexesWetSurfaces) ;
    
    nsurf = length(IndexesWetSurfaces) ;  % Number of wet surfaces
    
    
    IndexElemSelectVECT = cell2mat(IndexElemSelect') ;
    
    CNb =  MESH.CNb(IndexElemSelectVECT,:)  ;   % These are the connectivity matrix
    % of the nodes belonging to the
    % wet surfaces
    
    MaterialLibraryBND  = (nsurf+1)*ones(size(CNb,1),1) ;
    iini = 1;
    for isurfLOC  = 1:nsurf
        nelemsLOC = length(IndexElemSelect{isurfLOC}) ;
        ifin = iini + nelemsLOC-1;
        MaterialLibraryBND(iini:ifin) = IndexesWetSurfaces(isurfLOC) ;
        iini = ifin+1 ;
    end
    
end

MESH.HYDRO = [] ; 
MESH.HYDRO.CNb = CNb ; 


% 1. Assembly operator \LLbB{}{}  (Gauss point level )
% ---------------------------------


nelemB = size(CNb,1) ; % Number elementes boundary
nnodeE = size(CNb,2) ; % Number of nodes per element (boundary)

DATA.MESH.HYDRO = [] ; 
DATA.MESH.HYDRO.nelemB = nelemB ; 


TypeIntegrand = 'RHS';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(MESH.TypeElementB,nnodeE,TypeIntegrand) ;

ngaus = size(posgp,2) ;   % Number of Gauss points.
nnode = size(MESH.COOR,1) ; % Number of nodes, for all the discretization
ndim = 3;
Lbool = Lbool_vectorized(CNb,nnode,ngaus,ndim)  ;
%Lbool = Lbool_nonvectorized(CNb,nnode,ngaus,ndim)  ;

DATA.MESH.HYDRO.nelemB = nelemB ; % Number of elements wet surface 
DATA.MESH.HYDRO.ngausT = nelemB*ngaus ; % Number of elements wet surface 
DATA.MESH.HYDRO.ngaus = ngaus ;  
DATA.MESH.HYDRO.nnodeE = nnodeE ; 
% ------------------------------------------------------
% 2. MAtrix of shape functions for ALL Gauss points (concatenation, and multiplied by the weights)
% -------------------------------------------------
wSTft =repmat(weig',nelemB,1) ;      ; % Vector of weights for all Gauss points 


Nescl_w = bsxfun(@times,shapef,weig') ;  % Scalar functions multiplied by the weights

% Each entry multiplied by the identiy  N1*IDENT, N2*IDENT2 ....
Nb_xi_w = zeros(ndim*size(Nescl_w,1),ndim*size(Nescl_w,2)) ;
for igaus = 1:size(Nescl_w,1)
    INI = (igaus-1)*ndim+1; FIN = igaus*ndim ;
    Nb_xi_w(INI:FIN,:) = StransfN(Nescl_w(igaus,:),ndim) ;
end

% Since we are assuming that all boundary elements are of the same family,
% this matrix consists in NGAUS_all tiling copies of shapef
% Strictly speaking, this would   not be necessary  in this
% case, as we are assuming that all elements pertains to the same family (we do it for the sake
% of tenerality)
NbST_w_elem = repmat(Nb_xi_w, nelemB,1) ;



% ---------------------------------
% 3. Convert NbST_w_elem  into a diagonal matrix
%
% Now we have NbST_w_elem = [N_1; N_2 ... N_M]  (M gauss points), and we
% wish to turn it into a block diagonal matrix
% We have to adapt /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/ConvertBlockDiag.m% to cope
%   *** /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/ConvertBlockDiag_general.m
[NbST_diag_w,irowsNDIM,icolsNDIM ]= ConvertBlockDiag_general(NbST_w_elem,ndim) ;

NbST_w =  NbST_diag_w*Lbool ;


% 4.  Matrix BstB_allgauss, such that   tau_i = diag(BstB_allgauss{i})*Lbool*COOR'(:) gives the
% tangent vectors in the i-direction
% \BstB{i}{e}
% --------------------------------------------------------------------------
% For a given element
% -------------------
BstB_e = cell(2,1) ;
BstB_e{1} = zeros(ndim*size(Nescl_w,1),ndim*size(Nescl_w,2)) ;
for igaus = 1:size(Nescl_w,1)
    INI = (igaus-1)*ndim+1; FIN = igaus*ndim ;
    for  idim  = 1:2
        BstB_e{idim}(INI:FIN,:) = StransfN(squeeze(dershapef(idim,:,igaus)),ndim) ;
    end
end
% Since all elements are of the same family
BstB_allgauss = cell(2,1) ;
for  idim  = 1:2
    BstB_allgauss{idim}  = repmat(BstB_e{idim}, nelemB,1) ;
end

% -------------------------------------------------------------------------
% 5. Initial tangent vectors
% -------------------------
tTANGiniST = cell(1,2) ;
X = MESH.COOR' ;
for idim = 1:2
    tTANGiniST{idim} = ConvertBlockDiag_general(BstB_allgauss{idim},ndim,irowsNDIM,icolsNDIM)*Lbool*X(:) ;
end
% See assessment in FLOAT_IMPLE.pdf, kw:89


% 6. Initial normals
% ------------------
tau1 = reshape(tTANGiniST{1},ndim,[]) ;
tau2 = reshape(tTANGiniST{2},ndim,[]) ;
mINI =  cross(tau1,tau2) ;
mNORMiniST = mINI(:);

nmINI  = sqrt(sum(mINI.^2,1)) ;
%unitNORMALS = bsxfun(@times,mINI',1./nmINI')';


% 7. Component   of Nst in the gravity direction
% \NbSTgrav
%------------------- 
if DATA.vGRAVITY(2) == 0   
    error('GRAVITY SHOULD BE IN THE X2 DIRECTION')
end 
Nb_y = zeros(size(shapef,1),ndim*size(shapef,2)) ;
e2 = [0,1,0]' ; 
for igaus = 1:size(shapef,1)
   
    Nloc= StransfN(shapef(igaus,:),ndim) ;
    Nb_y(igaus,:) = e2'*Nloc ; 
   
end

NbSTgrav = repmat(Nb_y, nelemB,1) ;

% Initial x2 coordinate 
% \yPRESSiniST

[DDD,irows1,icols1] = ConvertBlockDiag_general(NbSTgrav,1) ; 

[yPRESSiniST]  = DDD*Lbool*X(:) ;


% 8. This is a little bit redundant (but we do it for the sake of clarity)
% Next we compute the shape functions of all boundary elements 

Nescl = shapef ;  % Scalar functions multiplied by the weights

% Each entry multiplied by the identiy  N1*IDENT, N2*IDENT2 ....
Nb_xi  = zeros(ndim*size(Nescl,1),ndim*size(Nescl,2)) ;
for igaus = 1:size(Nescl,1)
    INI = (igaus-1)*ndim+1; FIN = igaus*ndim ;
    Nb_xi(INI:FIN,:) = StransfN(Nescl(igaus,:),ndim) ;
end
 
NbST = repmat(Nb_xi, nelemB,1) ;


HYDRO.tTANGiniST = tTANGiniST ;   % Tangent vectors at each Gauss point (not unitary)
HYDRO.mNORMiniST = mNORMiniST ;   % Normal vector at each GAuss point (not unitary, it includes the Jacobian of the Gauss point )

HYDRO.BstB = BstB_allgauss;  % Derivate shape function at all Boundary Gauss points   

HYDRO.yPRESSiniST = yPRESSiniST;  % Initial coordinate X2 
HYDRO.Lbool = Lbool;  % Boolean operator 
HYDRO.NbSTgrav = NbSTgrav;  %  MAtrix such that NbSTgrav*Lbool*d gives the increment of x2
HYDRO.NbST_w = NbST_w;  %  MAtrix such that Fpr = -NbST_w^t*t_pr   
HYDRO.NbST = NbST;  %  
HYDRO.wST = wSTft ; 
HYDRO.JacobianWeights = nmINI ; % Jacobian at the Gauss point 
HYDRO.irows1 = irows1 ; 
HYDRO.icols1 = icols1 ; 


HYDRO.irowsNDIM = irowsNDIM ; 
HYDRO.icolsNDIM = icolsNDIM ; 

DATA.FOLLOWER_LOADS.HYDROSTATIC = DefaultField(DATA.FOLLOWER_LOADS.HYDROSTATIC,'WET_DAMPING_COEFFICIENT',0) ; %  = 1 ; 


    DATA.TMP.NO_COMPUTE_STIFFNESS_HYDRO = 0 ; 




%DRAW_NORMALS = 0;

% if DRAW_NORMALS == 1
%     
%     NameFile_msh = [DATA.NAMEgidFOLDER,filesep,'NORMALS','.msh'] ;
%     NameFile_res = [DATA.NAMEgidFOLDER,filesep,'NORMALS','.res'] ;
%     
%     ntau1  = sqrt(sum(tau1.^2,1)) ;
%     tau1 = bsxfun(@times,tau1',1./ntau1')';
%     
%     ntau2  = sqrt(sum(tau2.^2,1)) ;
%     tau2 = bsxfun(@times,tau2',1./ntau2')';
%     

%     
%     RotMatrixGlo = zeros()
%     
%     GidMesh2DFE_multi(NameFile_msh,MESH.COOR,{CNb},'',{MaterialLibraryBND},{MESH.TypeElementB},{'BND'});
%     NormalsPLOT(DATAIN,{MESH.TypeElementB},{'BND'},NameFile_res,...
%         tau1,tau2,mINI) ;
% end



