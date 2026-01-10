function [weig,posgp,shapef,dershapef,COOR] = Hexahedra27NInPoints(TypeIntegrand) 
%%
% =========================================================================
% Hexahedra27NInPoints — Gauss rule & shape functions for 27-node (triquadratic) hexahedral element
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, Q2 (triquadratic) shape functions, and their
%   derivatives w.r.t. parent coordinates (ξ, η, ζ) for a 27-node isoparametric
%   hexahedron. By default, uses a 3×3×3 Gauss–Legendre rule; a custom rule may
%   be injected by the caller.
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef, COOR] = Hexahedra27NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS'  (both use the default 3×3×3 Gauss–Legendre rule)
%     • cell  : {posgp, weig} to override the default rule
%               - posgp : 3×ngaus parent coordinates [ξ; η; ζ]
%               - weig  : 1×ngaus Gauss weights
%   (If nargin==0, the function loads 'tmp2.mat' for debugging/regression tests.)
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 3×ngaus Gauss points (ξ, η, ζ) in the parent domain
%   shapef      : ngaus×27 matrix of Q2 shape functions N(ξ,η,ζ)
%   dershapef   : 3×27×ngaus array with [∂N/∂ξ; ∂N/∂η; ∂N/∂ζ] at each Gauss point
%   COOR        : 27×3 array of nodal parent coordinates for the 27-node element
%
% DEFAULT QUADRATURE (3×3×3 Gauss–Legendre)
%   • Points  : ξ,η,ζ ∈ {−√(3/5), 0, √(3/5)} (tensor product grid → 27 points)
%   • Weights : {5/9, 8/9, 5/9} (tensor product to obtain 27 weights)
%   The routine COOR_Hexag27(·) provides both nodal coordinates (axy=bz=1)
%   and the Gauss grid (axy=bz=√(3/5)).
%
% IMPLEMENTATION NOTES
%   • Node numbering follows GiD: 8 corners, 12 midsides, 6 face centers, 1 center.
%     See: https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type
%   • Shape functions are constructed via a tensor-product quadratic polynomial basis
%     with orders n = [2 2 2]. The code forms Vandermonde-like matrices:
%         Pnodes = CoordinateMatrixPolynomial(COOR, n)
%         COEFFS = inv(Pnodes)
%         Pgauss, Pder = CoordinateMatrixPolynomial(posgp, n)
%     Then:
%         shapef          = Pgauss * COEFFS               % (ngaus×27)
%#         dershapef(idim,:,:)= (Pder{idim} * COEFFS)'     % idim = 1..3 → ξ,η,ζ
%   • dershapef is given in parent coordinates; mapping to physical (x,y,z)
%     requires the Jacobian of the isoparametric transformation at each Gauss point.
%   • Custom quadrature: pass {posgp, weig} with posgp as 3×ngaus and weig as 1×ngaus.
%
% USAGE
%   [w, xg, N, dN, Xn] = Hexahedra27NInPoints('K');             % stiffness/stress assembly
%   [w, xg, N, dN, Xn] = Hexahedra27NInPoints('RHS');           % mass/loads assembly
%   [w, xg, N, dN, Xn] = Hexahedra27NInPoints({Xc, Wc});        % custom quadrature
%
% REFERENCES (GiD conventions)
%   • Element type & numbering:
%     https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type
%   • (Related) Gauss points in GiD post format:
%     https://www.gidhome.com/documents/customizationmanual/POSTPROCESS%20DATA%20FILES/Results%20format:%20ModelName.post.res/Gauss%20Points
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

% This function returns, for each 27-node quadrilateral elements,
% and using a 3x3x3 Gauss rule (ngaus=0), or a Gauss rule given by the user
% NUMBERING CONVENTION: https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type
% weig = Vector of Gauss weights (1xngaus)
% posgp: Position of Gauss points  (ndim x ngaus)
% shapef: Array of shape functions (ngaus x nnodeE)
% dershape: Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)
% Joaquin A. Hernandez,  10th-April-2020.  (28th day of COVID-19 quarantine )
% Revised: 7th-January-2022 (see )
%-----------------------------------------------------------------------------

%[N,B] = Hexahedra27points_symbolic ;

%addpath(genpath('/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/SVDlibrary')) ;
if nargin == 0
    load('tmp2.mat')
end

nnodeE= 27;
ndim = 3;
% COORDINATES OF THE NODES
% ------------------------
axy = 1;
bz = 1;
COOR = COOR_Hexag27(axy,bz) ;
if ~iscell(TypeIntegrand)
    % COORDINATES OF THE GAUSS POINTS
    axy = sqrt(3/5) ;
    bz = axy ;
    [posgp,weig ]= COOR_Hexag27(axy,bz) ;
else
    weig = TypeIntegrand{2} ;
    posgp = TypeIntegrand{1}' ;
end

% Now we wish to write  the shape function of node I as
% NI = (x^0*y^0*z^0)*a_1 +(x^1*y^0*z^0)*a_2 + .... (x^2*y^2*z^2)*a_27
%
n = [2 2 2] ;
[Pnodes ]= CoordinateMatrixPolynomial(COOR,n) ;
COEFFSpol = inv(Pnodes) ;

% Shape functions at the Gauss Points
[Pgauss,Pder ]= CoordinateMatrixPolynomial(posgp,n) ;
shapef = Pgauss*COEFFSpol ;    % ngaus x nnodeE
ndim = 3;  nnodeE = 27 ; ngaus = length(weig) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;

for idim = 1:ndim
    
    dershapef(idim,:,:) =  (Pder{idim}*COEFFSpol)' ;
    
end

posgp = posgp' ;
weig = weig' ;


end

