function [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionDIRECTFIT2023(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)
%--------------------------------------------------------------------------
% function [PHIk_y, dPHIk_y, POLYINFO] = EvaluateBasisFunctionDIRECTFIT2023( ...
%                    xNEW, DATA, VAR_SMOOTH_FE, POLYINFO)
%
% PURPOSE:
%   Evaluates local basis functions and their derivatives at arbitrary
%   spatial points (`xNEW`) by identifying the element that contains each
%   point and constructing a **direct-fit polynomial interpolant**.
%   Supports applications such as empirical interpolation, smoothing,
%   and post-processing.
%
% INPUTS:
% -------
%   - xNEW          : [npoints × ndim] array of coordinates where basis
%                     functions are to be evaluated.
%   - DATA          : structure controlling behavior of evaluation.
%                     · OnlyCheckIfIsInside → (bool) if 1, only checks if
%                       points are inside mesh (no basis function eval)
%   - VAR_SMOOTH_FE : structure with fields:
%                     · CN                  → connectivity matrix
%                     · COOR                → nodal coordinates
%                     · ORDER_POLYNOMIALS   → order of local polynomials
%                     · BasisIntegrand      → coefficients for integrands
%                     · ExactIntegral       → reference values
%   - POLYINFO      : auxiliary structure storing:
%                     · Previous triangulations
%                     · Polynomial coefficients
%                     · Scaling for local variables
%
% OUTPUTS:
% --------
%   - PHIk_y   : [npoints × nfun] matrix with values of basis functions at
%                each query point.
%   - dPHIk_y  : 1×ndim cell, each cell is [npoints × nfun] and contains
%                directional derivatives of basis functions.
%   - POLYINFO : updated structure with:
%                · setElements: element index for each query point
%                · COEFFSpolynomial: fitted polynomial coefficients
%                · ListPointsOutside: points not found in the mesh
%
% METHOD:
% -------
% For each query point:
%   1. Identify closest node and candidate elements.
%   2. Use `WhichElementInside2023` to locate the containing element.
%   3. If inside:
%       a. Retrieve or compute polynomial coefficients via
%          `CoeffsPolyShapeFunctions2023`.
%       b. Evaluate basis function and derivatives using
%          `CoordinateMatrixPolynomial`.
%   4. If outside:
%       - Log in POLYINFO.ListPointsOutside.
%
% EXAMPLE:
% --------
%   [phi, dphi, poly] = EvaluateBasisFunctionDIRECTFIT2023(Xeval,DATAin,VAR_FE,POLY0);
%
% REMARKS:
% --------
%   - Polynomial fitting is local and reused via `POLYINFO`.
%   - Coordinates are normalized element-wise for numerical stability.
%   - Points outside the domain return empty entries.
%
% DEPENDENCIES:
% -------------
%   · `WhichElementInside2023`
%   · `FindClosestNodes`
%   · `CoeffsPolyShapeFunctions2023`
%   · `CoordinateMatrixPolynomial`
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, March 2023
%--------------------------------------------------------------------------



if nargin == 0
    load('tmp2.mat')
end

DATA = DefaultField(DATA,'OnlyCheckIfIsInside',0) ;

npoints = size(xNEW,1) ; % Number of points at which to evaluate the function
% Number of functions to be integrated
nfun = length(VAR_SMOOTH_FE.ExactIntegral) ;
% Initialization outputs
PHIk_y =zeros(npoints,nfun) ;  % Value of the integrand function at xNEW.
ndim =  size(VAR_SMOOTH_FE.COOR,2) ; % Number of spatial dimensions
dPHIk_y = cell(1,ndim) ;  % Value of the gradient of the integrand functions at xNEW
for idim = 1:ndim
    dPHIk_y{idim} = zeros(npoints,nfun) ;   % Initialization with zeros
end
nelem = size(VAR_SMOOTH_FE.CN,1)  ;  % Total number of FE elements
nnodeE = size(VAR_SMOOTH_FE.CN,2)  ;  % Number of nodes per element

% Coefficients of the interpolation polynomial for each element
% At the beginning, they are empty. As the search algorithm proceeds,
% coefficients are calculated and stored in this cell array.
POLYINFO = DefaultField(POLYINFO,'COEFFSpolynomial',cell(nelem,1)) ;%
% Find the nodes of the FE discretization which are closer to xNEW
[INDEXES_NEAR,IND_POLYG] = FindClosestNodes(xNEW,VAR_SMOOTH_FE) ;
ELEMENTS_CONTAINING_xNEW = zeros((npoints),1) ; % Indices of the elements containing the points
%
%INDfin = 0 ;
inew = 1;
POLYINFO = DefaultField(POLYINFO,'SCALING_VARIABLES',[]) ;
POLYINFO.SCALING_VARIABLES = DefaultField(POLYINFO.SCALING_VARIABLES,'LENGTH',zeros(nelem,ndim)) ;
POLYINFO.SCALING_VARIABLES = DefaultField(POLYINFO.SCALING_VARIABLES,'coorREF',zeros(nelem,ndim)) ;

%SCALING_ACTIVE = 1;
%BdomRED_interp = cell(npoints,1) ;  % Interpolation of the reduced B-matrix
POLYINFO.ListPointsOutside = [] ; % Lists points outside

while inew <=npoints
    xLOC = xNEW(inew,:); % Coordinate of the point at which the function is to be evaluated
    INDnear = INDEXES_NEAR(inew) ;  % INDEX Nearest FE mesh  node to xLOC
    % Searching for the element containing xLOC --> elemCONTAINER
    [elemCONTAINER,POLYINFO ]= WhichElementInside2023(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG,POLYINFO,inew) ;
    
%     [elemCONTAINER,POLYINFO ]= WhichElementInside(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG,POLYINFO,inew) ;
    
    %  elemCONTAINER = WhichElementInside(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG) ;
    if isempty(elemCONTAINER)
        % Point outside the domain
        PHIk_y = [];
        ELEMENTS_CONTAINING_xNEW(inew ) = 0 ;
        POLYINFO.ListPointsOutside = [POLYINFO.ListPointsOutside;inew] ;
%         if DATA.OnlyCheckIfIsInside == 0
%             disp('')
%         end
        
    else
        ELEMENTS_CONTAINING_xNEW(inew ) = elemCONTAINER ;  % Element containing the point
    end
    
    
    if DATA.OnlyCheckIfIsInside == 0   && ~isempty(elemCONTAINER)
     %     if DATA.OnlyCheckIfIsInside == 0  % Before march-17-2023
        % Determining the matrices for constructing shape functions
        NODESloc = small2large(elemCONTAINER,VAR_SMOOTH_FE.ngausE) ;
        [COEFFSpol,LelemSCALING,coorREF,POLYINFO] = CoeffsPolyShapeFunctions2023(elemCONTAINER,NODESloc,VAR_SMOOTH_FE,POLYINFO)  ;
        
        % Evaluating the shape functions (and derivatives) at  point xLOC
        % ---------------------------------------------------------------------
        xLOC  = (xLOC-coorREF)./LelemSCALING   ;  % Transformed coordinate of the point under study
        factorDERIVATIVE = 1./LelemSCALING ;    % This for computing the derivatives
        % -------------------------------------------------------------------------------------------
        % Computatio  of the shape functions (recall that the "nodes" here of such functions are the Gauss points)
        [Pevaluate, PevaluateDER]= CoordinateMatrixPolynomial(xLOC,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
        N = Pevaluate*COEFFSpol ;  % Shape functions at the given points COORevaluate
        PHIk_y(inew,:)  =   N*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
        for idim = 1:length(PevaluateDER)
            dN = PevaluateDER{idim}*COEFFSpol ;
            dPHIk_y{idim}(inew,:) = factorDERIVATIVE(idim)*dN*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
        end        
    end
    
    inew = inew +1 ;
end

POLYINFO.setElements = ELEMENTS_CONTAINING_xNEW ; %(inew ) = elemCONTAINER ;  % Element containing the point

%POLYINFO.ELEMENTS_CONTAINING_xNEW = ELEMENTS_CONTAINING_xNEW ; %(inew ) = elemCONTAINER ;  % Element containing the point
%POLYINFO.BdomRED_interp = BdomRED_interp ;
