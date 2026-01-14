function D_PdownsRBlROT_e = D_PdownsRBlROT_COROT(EIFEoper_all,lambdaLEN,DOFsROT,DOFsTR)
%--------------------------------------------------------------------------
%  D_PdownsRBlROT_COROT
%
%  Computes the element-wise **tangent downscaling operator** relating 
%  variations of the **rotation vector** (in local coordinates) to variations
%  of the **coarse-scale displacement vector**, under the **co-rotational EIFEM**
%  formulation.
%
%  OUTPUT:
%    - D_PdownsRBlROT_e : Matrix used to express δθ (variation of rotation vector) as:
%
%           δθ = P̂γr Qᵀ δd_c     (see Eq. (398)):contentReference[oaicite:0]{index=0}
%
%  THEORY:
%    This matrix corresponds to the rotational component of the small-strain
%    downscaling operator for rigid-body displacements (last rows of P̂), normalized
%    by the element's scale parameter λ (lambdaLEN), i.e.:
%
%           D_PdownsRBlROT_e = (1 / λ) * PdownsRB_rot    (Eq. 582):contentReference[oaicite:1]{index=1}
%
%    where PdownsRB is computed during the **training phase** and stored in:
%       > EIFEoper_all.OPER.PdownsRB
%
%    The required block is extracted via the indices `DOFsROT` (rotation) and
%    `DOFsTR` (translation). The function also checks that translational and rotational
%    modes are properly decoupled (i.e., origin at interface centroid):
%
%           norm(Grb_TR) / norm(Grb_ROT) ≪ 1
%
%    Failure of this check indicates an inconsistent reference system, and stops execution.
%
%  APPLICATION:
%    This matrix is critical for:
%     - Computing geometric tangent stiffness (e.g., via `Ssp(F_int) * P̂γr`)
%     - Relating δQ with δd_c (see Sections 12.6 and 12.7)
%     - Assembly of `Kc_geo` in both 2D and 3D:contentReference[oaicite:2]{index=2}
%
%  INPUTS:
%    - EIFEoper_all : structure with stored Grb and PdownsRB from training
%    - lambdaLEN    : element-wise scale factor (computed in parent domain mapping)
%    - DOFsROT      : indices of rotational DOFs
%    - DOFsTR       : indices of translational DOFs
%
%  REMARK:
%    This operator ensures **objectivity** of the tangent formulation by mapping variations
%    from physical displacements to a local rotational frame.
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 7-Feb-2025
%    Comments by ChatGPT4, 13-May-2025
%
%  Related expressions:
%    - Eq. (255), (256), (258), (582) in *EIFEM_largeROTfinal.pdf*
%    - P̂γr = [P̂γrb P̂γr0], used in geometric stiffness: K_geo = Ssp(F_int) P̂γr Qᵀ
%     
%--------------------------------------------------------------------------





% Tangent downscaling matrix relating variations of the rotation angle (in local coordinates)
%with variations of  the vector of coarse-scale displacements
%($\DiagC{\PdownsRBlROT}$,
% JAHO, 7-FEB-2025, 8:01, Balmes 185, Barcelona

% In  EIFE_operatorsBUB.m, we computed
% EIFEoper.OPER.Grb = Grb;
% EIFEoper.OPER.PdownsRB = PdownsRB;
% First we need to be sure that Grb is "decoupled"
GrbTRrot = EIFEoper_all.OPER.Grb(DOFsTR,DOFsROT) ;
GrbROTrot = EIFEoper_all.OPER.Grb(DOFsROT,DOFsROT) ;
coupligROTtr = norm(GrbTRrot)/norm(GrbROTrot) ;
if coupligROTtr > 1e-10
    disp('There is coupling between rotation and translations')
    disp('Check that, in the offline stage, the origin of the local coordinate system  ')
    disp('coincides with the centroid of the boundary interfaces')
    error('')
end
%\PdownsRBlROTe{e} = \dfrac{1}{\lambdaLENe{e}}  \PdownsRBrotE{e}
% \PdownsRBrotE{e}$ is the rotational component (last $\nRB$ rows) of the small strains  downscaling operator
% for rigid body displacements
D_PdownsRBlROT_e = sparse(EIFEoper_all.OPER.PdownsRB(DOFsROT,:)/lambdaLEN);