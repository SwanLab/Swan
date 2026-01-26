function [D_AspMAT_e,DOFsROT,DOFsTR] =  SpinOperator(ndim,nnodeCb,NumberBubbleDOFS_perELEMENT_e)
%--------------------------------------------------------------------------
%  SpinOperator
%
%  Constructs the **spin operator matrix** used in the co-rotational EIFEM 
%  framework to represent infinitesimal rotations in 2D and 3D problems.
%
%  The output matrix `D_AspMAT_e` plays a central role in the geometric 
%  stiffness component of the tangent matrix and in the consistent variation
%  of the rotation matrices, especially in small-strain/large-rotation regimes.
%
%  Mathematical background:
%    For a node at position X, the infinitesimal rotation is expressed using the
%    skew-symmetric **spin tensor**, defined in 3D as:
%
%      spin(u) = [  0  -u3  u2
%                  u3   0  -u1
%                 -u2  u1   0 ]   :contentReference[oaicite:0]{index=0}
%
%    For 2D, this reduces to:
%      ssp = [ 0  -1
%              1   0 ]
%
%    When applied across all nodes, this yields a block-diagonal matrix with 
%    spin blocks on the diagonal (one per boundary node). Bubble DOFs are 
%    not affected by the rotation and are assigned identity blocks.
%
%  OUTPUTS:
%    - D_AspMAT_e : 
%         > 2D: a single block-diagonal matrix (sparse) combining spin operators
%         > 3D: a cell array {D_AspMAT_x, D_AspMAT_y, D_AspMAT_z}
%               where each entry corresponds to the spin matrix for one rotation axis
%
%    - DOFsROT : Global indices of rotational DOFs (used in downscaling operator)
%    - DOFsTR  : Global indices of translational DOFs
%
%  Role in the formulation:
%    - `D_AspMAT_e` is used to compute geometric stiffness contributions:
%         δQ ≈ spin(δar) Q   ⇒ geometric variations involve `spin()` matrices:contentReference[oaicite:1]{index=1}
%    - In Section 12.6.3, this spin matrix appears in the assembly of `Ssp(F_int)`
%      and the geometric stiffness terms of the form:
%         Dn(δQ) F_int = Ssp(F_int) P̂γr Qᵀ δd_c  :contentReference[oaicite:2]{index=2}
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 24-Feb-2024
%    Comments by ChatGPT4, 13-May-2025
%  Related implementations:
%    - Used in `B_N_matricesEIFEbubCOROT_LRss.m`
%    - Appears in `D_PdownsRBlROT_COROT.m` for computing downscaling operators
%    - Theoretical formulation in Sections 12.6–12.8 and Eq.(590) of *EIFEM_largeROTfinal.pdf*
%
%--------------------------------------------------------------------------




% EIFE METHOD, DETERMINATION OF SPIN OPERATOR
% CO-ROTATIONAL APPROACH,
% /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/06_COROT_3D.mlx
% JAHO: 24-Feb-2O24,  Monday, 7:39 AM,  BALMES 185, BARCELONA
if nargin == 0
    load('tmp1.mat')
end
if  ndim == 2
    % --------------------------
    % SPIN OPERATOR
    % --------------------------
    %   \AspMATe{e} \defeq  \matcdos{\sSPmatE{e}}{\zero}{\zero}{\zero}
    %    \sSPmatE{e} \defeq  \diagOL{\sSPmat}{\sSPmat}{\cdots}{\sSPmat}
    %    \sSPmat = \matcdos{0}{-1}{1}{0}
    sSPmatE = cell(nnodeCb,1) ;
    sSPmatE(:) =  {sparse([0 -1;
        1  0])} ;
    sSPmatE = blkdiag(sSPmatE{:}) ;
    AspMATe = cell(2,1) ;
    AspMATe{1} = sSPmatE ;
    AspMATe{2} = sparse(NumberBubbleDOFS_perELEMENT_e,NumberBubbleDOFS_perELEMENT_e) ;
    D_AspMAT_e = blkdiag(AspMATe{:}) ;
    
    DOFsROT = 3;
    DOFsTR = 1:2 ;
else
    
    ndim = 3; 
    sSPmatE = cell(1,ndim) ;
    %   \sSPmatTD{1} = \matctres{0}{0}{0}
    %   {0}{0}{-1}
    %   {0}{1}{0},
    idim = 1;
    sSPmatE{idim} = cell(nnodeCb,1) ;
    sSPmatE{idim}(:) =  {sparse([0 0  0;
        0  0   -1 ;
        0  1    0])} ;
    %\sSPmatTD{2} = \matctres{0}{0}{1}
    %   {0}{0}{0}
    %   {-1}{0}{0},
    idim = 2;
    sSPmatE{idim} = cell(nnodeCb,1) ;
    sSPmatE{idim}(:) =  {sparse([0 0  1;
        0  0   0 ;
        -1  0    0])} ;
    %   \sSPmatTD{3} = \matctres{0}{-1}{0}
    %   {1}{0}{0}
    %   {0}{0}{0}
    idim = 3;
     sSPmatE{idim} = cell(nnodeCb,1) ;
     sSPmatE{idim}(:)=  {sparse([0 -1  0;
        1  0   0 ;
        0  0    0])} ;
    
    %  \AspMATtdE{e}{i} \defeq  \matcdos{\sSPmatTDe{e}{i}}{\zero}{\zero}{\zero}, \hspace{0.5cm}  i = 1,2,3
  %AspMATe = cell(1,ndim) ; 
   D_AspMAT_e = cell(1,ndim) ; 
   for idim = 1:ndim
       AspMATe =  cell(2,1) ; 
         sSPmatE{idim}  = blkdiag(sSPmatE{idim}{:}) ;
       AspMATe{1} = sSPmatE{idim} ; 
       AspMATe{2} = sparse(NumberBubbleDOFS_perELEMENT_e,NumberBubbleDOFS_perELEMENT_e) ;
       D_AspMAT_e{idim} = blkdiag(AspMATe{:}) ;
   end  
   
    
    
    DOFsROT = 4:6;
    DOFsTR = 1:3 ;
end