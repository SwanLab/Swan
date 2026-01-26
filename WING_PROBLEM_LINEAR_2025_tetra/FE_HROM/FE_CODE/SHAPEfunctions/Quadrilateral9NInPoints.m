function [weig,posgp,shapef,dershapef] = Quadrilateral9NInPoints(TypeIntegrand)
%%
% =========================================================================
% Quadrilateral9NInPoints — Gauss rule & shape functions for 9-node (quadratic) quadrilateral
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, quadratic (Q2) shape functions, and their
%   derivatives w.r.t. parent coordinates (ξ, η) for a 9-node isoparametric
%   quadrilateral. By default, uses a 3×3 Gauss–Legendre rule, suitable for
%   stiffness and mass/load integrals in isoparametric Q2 elements.
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Quadrilateral9NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS'  (both default to 3×3 Gauss–Legendre)
%     • cell  : {posgp, weig} to override the default rule
%               - posgp : 2×ngaus parent coordinates [ξ; η]
%               - weig  : 1×ngaus Gauss weights
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 2×ngaus Gauss points (ξ, η in parent domain)
%   shapef      : ngaus×9 matrix of Q2 shape functions N(ξ,η)
%   dershapef   : 2×9×ngaus ∂N/∂ξ, ∂N/∂η at each Gauss point
%
% DEFAULT QUADRATURE (3×3 Gauss–Legendre)
%   • Points  : ξ,η ∈ {−√(3/5), 0, √(3/5)} (tensor product grid)
%   • Weights : {5/9, 8/9, 5/9} (tensor product)
%   The code internally reorders the 9 points with CONV to follow the
%   conventional GiD order.
%
% IMPLEMENTATION NOTES
%   • Node numbering (GiD convention): corners 1–4 (from (−1,−1) CCW),
%     midsides 5–8, center 9.
%   • Shape functions N(ξ,η) and derivatives B(ξ,η) are defined as local
%     function handles (N, B). B returns:
%         BeXi = [∂N/∂ξ; ∂N/∂η]  (size 2×9)
%   • dershapef is given in parent coordinates. Mapping to physical (x,y)
%     requires the Jacobian of the isoparametric mapping.
%   • Custom quadrature can be injected via TypeIntegrand = {posgp, weig}.
%
% USAGE
%   [w, xg, N, dN] = Quadrilateral9NInPoints('K');          % stiffness assembly
%   [w, xg, N, dN] = Quadrilateral9NInPoints('RHS');        % mass/loads assembly
%   [w, xg, N, dN] = Quadrilateral9NInPoints({Xc, Wc});     % custom rule
%
% REFERENCES (GiD format and conventions)
%   • Preprocess element types:
%     https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20typ
%   • Postprocess Gauss points:
%     https://www.gidhome.com/documents/customizationmanual/POSTPROCESS%20DATA%20FILES/Results%20format:%20ModelName.post.res/Gauss%20Points
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================


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



N = @(x,y) [x.*(x-1).*y.*(y-1)/4, x.*(x+1).*y.*(y-1)/4, ...
    x.*(x+1).*y.*(y+1)/4, x.*(x-1).*y.*(y+1)/4, ...
    (1-x.^2).*y.*(y-1)/2,  x.*(x+1).*(1-y.^2)/2,   ...
    (1-x.^2).*y.*(y+1)/2,  x.*(x-1).*(1-y.^2)/2,   ...
    (1-x.^2).*(1-y.^2)];



B = @(x,y) [(x-1/2).*y.*(y-1)/2,   (x+1/2).*y.*(y-1)/2, ...
    (x+1/2).*y.*(y+1)/2,   (x-1/2).*y.*(y+1)/2, ...
    -x.*y.*(y-1),          (x+1/2).*(1-y.^2),   ...
    -x.*y.*(y+1),          (x-1/2).*(1-y.^2),   ...
    -2*x.*(1-y.^2);
    x.*(x-1).*(y-1/2)/2,    x.*(x+1).*(y-1/2)/2, ...
    x.*(x+1).*(y+1/2)/2,    x.*(x-1).*(y+1/2)/2, ...
    (1-x.^2).*(y-1/2),       x.*(x+1).*(-y),   ...
    (1-x.^2).*(y+1/2),       x.*(x-1).*(-y),   ...
    (1-x.^2).*(-2*y)];


shapef = zeros(length(weig),9);
for i=1:length(weig)
    shapef(i,:) = N(posgp(1,i),posgp(2,i));
end

dershapef = zeros(2,9,length(weig));

for j=1:length(weig)
    dershapef(:,:,j) = B(posgp(1,j),posgp(2,j));
end
end