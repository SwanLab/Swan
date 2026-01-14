function [COOR,W ]= COOR_Hexag27(a,b)
%%
% =========================================================================
% COOR_Hexag27 — Node and Gauss point coordinates for 27-node (triquadratic) hexahedral element
% =========================================================================
% PURPOSE
%   Generate nodal or Gauss-point coordinates for a 27-node (Q2) hexahedral element
%   in parent coordinates (ξ, η, ζ), following the GiD numbering convention.
%   When called with a = b = 1 → nodal coordinates.
%   When called with a = b = √(3/5) → 3×3×3 Gauss–Legendre integration grid.
%
% SIGNATURE
%   [COOR, W] = COOR_Hexag27(a, b)
%
% INPUT
%   a, b : Scalar parameters controlling the coordinates along ξ, η, ζ axes.
%          Typically:
%             a = b = 1           → nodal positions (−1, 0, +1)
%             a = b = √(3/5)      → Gauss–Legendre positions for 3×3×3 rule
%
% OUTPUT
%   COOR : 27×3 matrix with coordinates [ξ, η, ζ] of the 27 nodes/points
%   W    : 27×1 vector with corresponding Gauss weights
%
% NUMBERING CONVENTION
%   The numbering follows the GiD reference standard:
%     • 1–4  : corner nodes (z = −1 plane)
%     • 5–8  : corner nodes (z = +1 plane)
%     • 9–12 : midside edges (z = −1)
%     • 13–16: midside edges (z = 0)
%     • 17–20: midside edges (z = +1)
%     • 21   : center (z = −1)
%     • 22–25: midside edges (z = 0)
%     • 26   : center (z = +1)
%     • 27   : element centroid (z = 0)
%
% FORMULATION DETAILS
%   • The function defines the coordinates layer by layer in z:
%       z = −b → bottom plane
%       z =  0 → mid plane
%       z = +b → top plane
%   • Each node is assigned its coordinate and local weight W(i),
%     using tensor-product Gauss–Legendre weights:
%         5/9, 8/9, 5/9  → along each dimension
%     so that:
%         W = wξ * wη * wζ
%   • When used for nodal coordinates (a=b=1), the weights W are irrelevant.
%
% USAGE
%   [Xn, Wn] = COOR_Hexag27(1,1);          % nodal coordinates (−1,0,+1)
%   [Xg, Wg] = COOR_Hexag27(sqrt(3/5),1);  % 3×3×3 Gauss–Legendre grid
%
% REFERENCES
%   • GiD element type (Hexahedra27):
%     https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type
%   • GiD Gauss points convention:
%     https://www.gidhome.com/documents/customizationmanual/POSTPROCESS%20DATA%20FILES/Results%20format:%20ModelName.post.res/Gauss%20Points
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

COOR = zeros(27,3) ;
W = zeros(27,1) ;
% NODES 1-2-3-4: Corners plane z = -1

z = -b;  wLOC =  (5/9)^3 ;
COOR(1,:) = [-a,-a,z] ;  W(1) = wLOC   ;
COOR(2,:) = [+a,-a,z] ;  W(2) = wLOC   ;
COOR(3,:) = [+a,+a,z] ;  W(3) =  wLOC  ;
COOR(4,:) = [-a,+a,z] ;  W(4) =  wLOC  ;
% NODES 5-6-7-8: Corners plane z = +1
z = +b;
COOR(5,:) = [-a,-a,z] ;  W(5) = wLOC ;
COOR(6,:) = [+a,-a,z] ;  W(6)  = wLOC ;
COOR(7,:) = [+a,+a,z] ;  W(7)  = wLOC ;
COOR(8,:) = [-a,+a,z] ;  W(8)  = wLOC ;
% NODES 9-10-11-12: Midside plane z = -1
z = -b;   wLOC =  (8/9)*(5/9)^2 ;
COOR(9,:) = [0,-a,z]  ;  W(9) = wLOC ;
COOR(10,:) = [+a,0,z] ;  W(10) = wLOC ;
COOR(11,:) = [0,+a,z] ;  W(11) = wLOC ;
COOR(12,:) = [-a,0,z] ;  W(12) = wLOC ;
% NODES 13-14-15-16 : Corners plane z = 0
z = 0;
COOR(13,:) = [-a,-a,z] ;  W(13) = wLOC ;
COOR(14,:) = [+a,-a,z] ;  W(14) = wLOC ;
COOR(15,:) = [+a,+a,z] ;  W(15) = wLOC ;
COOR(16,:) = [-a,+a,z] ;  W(16) = wLOC ;
% NODES 17-18-19-20: Midside plane z = +1
z =  +b;
COOR(17,:) = [0,-a,z]  ;   W(17) = wLOC ;
COOR(18,:) = [+a,0,z] ;   W(18) = wLOC ;
COOR(19,:) = [0,+a,z] ;   W(19) = wLOC ;
COOR(20,:) = [-a,0,z] ;  W(20) = wLOC ;
% NODES 21 - Center plane z = -1
z = -b ;
COOR(21,:) = [0,0,z] ;   W(21) = (5/9)*(8/9)^2 ;
% NODES 22-23-24-25: Midside plane z = +1
z =  0;    wLOC =  (5/9)*(8/9)^2  ;
COOR(22,:) = [0,-a,z]  ;   W(22) = wLOC ;
COOR(23,:) = [+a,0,z] ;    W(23) = wLOC ;
COOR(24,:) = [0,+a,z] ;    W(24) = wLOC ;
COOR(25,:) = [-a,0,z] ;   W(25) = wLOC ;
% NODES 26 - Center plane z = +1
z = +b ;
COOR(26,:) = [0,0,z] ;  W(26) = wLOC ;
% NODES 27 - Center plane z = 0
z = 0 ;
COOR(27,:) = [0,0,z] ;    W(27) = (8/9)^3 ;


end