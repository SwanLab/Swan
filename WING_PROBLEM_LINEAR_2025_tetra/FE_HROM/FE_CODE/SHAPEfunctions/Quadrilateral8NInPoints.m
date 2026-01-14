function [weig,posgp,shapef,dershapef] = Quadrilateral8NInPoints(TypeIntegrand)
%%
% =========================================================================
% Quadrilateral8NInPoints — Gauss rule & shape functions for 8-node (Q2-Serendipity) quadrilateral
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, serendipity-Q2 shape functions (8 nodes),
%   and their derivatives w.r.t. parent coordinates (ξ, η) for an isoparametric
%   quadrilateral. By default, a 3×3 Gauss–Legendre rule (tensor product) is used,
%   suitable for stiffness and mass/load integrals and robust for curved/distorted
%   meshes.
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Quadrilateral8NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS'  (both use the default 3×3 Gauss–Legendre rule)
%     • cell  : {posgp, weig} to override the default quadrature
%               - posgp : 2×ngaus parent coordinates [ξ; η]
%               - weig  : 1×ngaus Gauss weights
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 2×ngaus Gauss points (ξ, η in parent domain)
%   shapef      : ngaus×8 matrix of Q2-Serendipity shape functions N(ξ,η)
%   dershapef   : 2×8×ngaus array with [∂N/∂ξ; ∂N/∂η] at each Gauss point
%
% DEFAULT QUADRATURE (3×3 Gauss–Legendre)
%   • Points  : ξ,η ∈ {−√(3/5), 0, √(3/5)} on a tensor grid
%   • Weights : {5/9, 8/9, 5/9} (tensor product to obtain 9 weights)
%   The 9 points are reordered with CONV to follow the conventional GiD order.
%
% IMPLEMENTATION NOTES
%   • Node numbering (GiD convention): corners 1–4 (from (−1,−1) counter-clockwise),
%     midsides 5–8; this is an 8-node serendipity element (no center node).
%   • Shape functions and derivatives are obtained with the SEREN pipeline:
%         COOR_Quad8 → ShapeFunCoefficientsSEREN → ShapeFunctionFEseren
%     where ShapeFunctionFEseren returns N(ξ,η) and the derivative lists dN.
%   • dershapef is assembled from the returned derivative lists as:
%         dershapef(1,:,:) = ∂N/∂ξ,   dershapef(2,:,:) = ∂N/∂η
%   • Mapping to physical (x,y) requires the Jacobian of the isoparametric mapping.
%   • If nargin==0, the function loads 'tmp.mat' (handy for debugging/regression tests).
%
% USAGE
%   [w, xg, N, dN] = Quadrilateral8NInPoints('K');            % stiffness/stress assembly
%   [w, xg, N, dN] = Quadrilateral8NInPoints('RHS');          % mass/loads assembly
%   [w, xg, N, dN] = Quadrilateral8NInPoints({Xc, Wc});       % custom quadrature
%
% REFERENCES (GiD format and conventions)
%   • Preprocess element types:
%     https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20typ
%   • Postprocess Gauss points:
%     https://www.gidhome.com/documents/customizationmanual/POSTPROCESS%20DATA%20FILES/Results%20format:%20ModelName.post.res/Gauss%20Points
%   • See also:
%     /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/03_20nodeHEX.mlx
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

if nargin == 0
    load('tmp.mat')
end

if ~iscell(TypeIntegrand)
    % Standard 3-by-3 rule
    w1=5/9; w2=8/9; w3=5/9;
    weig= [w1*w1 w2*w1 w3*w1 w1*w2 w2*w2 w3*w2 w1*w3 w2*w3 w3*w3];
    p = sqrt(3/5);
    posgp=[-p   -p;    %   1
        0   -p;    %   5
        p   -p;    %   2
        -p    0;    %   8
        0    0;    %   9
        p    0;    %   6
        -p    p;    %   4
        0    p;    %   7
        p    p]';  %   3
    
    CONV = [1,5,2,8,9,6,4,7,3];
    
    [~,CONV] = sort(CONV) ;
    
    posgp = posgp(:,CONV) ;
    weig = weig(CONV) ;
else
    weig = TypeIntegrand{2}' ;
    posgp = TypeIntegrand{1} ; 
    
end


% SHAPE FUNCTIONS AND DERIVATIVE SHAPE FUNCTIONS 
% ...............--------------------------------.
ORDER_POLYNOMIALS = [2,2] ; 
[COORnodes]= COOR_Quad8 ; 
DATAshape = ShapeFunCoefficientsSEREN(COORnodes,ORDER_POLYNOMIALS) ;
DATAlocSHAPE.DATAshape  = DATAshape;
xLIM = [] ;
DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
[shapef,dN,~ ]=    ShapeFunctionFEseren(xLIM,posgp',DATAlocSHAPE) ;



dershapef = zeros(length(dN),size(shapef,2),length(weig));

 for idim=1:length(dN)
     dershapef(idim,:,:) =  dN{idim}' ; 
 end


% 
% 
% N = @(x,y) [x.*(x-1).*y.*(y-1)/4, x.*(x+1).*y.*(y-1)/4, ...
%     x.*(x+1).*y.*(y+1)/4, x.*(x-1).*y.*(y+1)/4, ...
%     (1-x.^2).*y.*(y-1)/2,  x.*(x+1).*(1-y.^2)/2,   ...
%     (1-x.^2).*y.*(y+1)/2,  x.*(x-1).*(1-y.^2)/2,   ...
%     (1-x.^2).*(1-y.^2)];
% 
% 
% 
% B = @(x,y) [(x-1/2).*y.*(y-1)/2,   (x+1/2).*y.*(y-1)/2, ...
%     (x+1/2).*y.*(y+1)/2,   (x-1/2).*y.*(y+1)/2, ...
%     -x.*y.*(y-1),          (x+1/2).*(1-y.^2),   ...
%     -x.*y.*(y+1),          (x-1/2).*(1-y.^2),   ...
%     -2*x.*(1-y.^2);
%     x.*(x-1).*(y-1/2)/2,    x.*(x+1).*(y-1/2)/2, ...
%     x.*(x+1).*(y+1/2)/2,    x.*(x-1).*(y+1/2)/2, ...
%     (1-x.^2).*(y-1/2),       x.*(x+1).*(-y),   ...
%     (1-x.^2).*(y+1/2),       x.*(x-1).*(-y),   ...
%     (1-x.^2).*(-2*y)];
% 
% 
% shapef = zeros(length(weig),9);
% for i=1:length(weig)
%     shapef(i,:) = N(posgp(1,i),posgp(2,i));
% end
% 
% dershapef = zeros(2,9,length(weig));
% 
% for j=1:length(weig)
%     dershapef(:,:,j) = B(posgp(1,j),posgp(2,j));
% end
% end