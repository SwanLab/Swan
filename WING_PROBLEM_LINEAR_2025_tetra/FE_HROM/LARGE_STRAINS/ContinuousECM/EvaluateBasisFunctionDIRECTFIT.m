function [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionDIRECTFIT(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)
% See
%  EvaluateBasisFunctionDIRECTFIT_aux.mlx


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
    elemCONTAINER = WhichElementInside(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG) ;     
    if isempty(elemCONTAINER)
       % Point outside the domain 
        PHIk_y = [];
        ELEMENTS_CONTAINING_xNEW(inew ) = 0 ;
        POLYINFO.ListPointsOutside = [POLYINFO.ListPointsOutside;inew] ; 
     
    else        
        ELEMENTS_CONTAINING_xNEW(inew ) = elemCONTAINER ;  % Element containing the point        
    end
    
    
    if DATA.OnlyCheckIfIsInside == 0
        
        % Determining the matrices for constructing shape functions
         NODESloc = small2large(elemCONTAINER,VAR_SMOOTH_FE.ngausE) ;
        [COEFFSpol,LelemSCALING,coorREF,POLYINFO] = CoeffsPolyShapeFunctions(elemCONTAINER,NODESloc,VAR_SMOOTH_FE,POLYINFO)  ; 
        
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


POLYINFO.ELEMENTS_CONTAINING_xNEW = ELEMENTS_CONTAINING_xNEW ; %(inew ) = elemCONTAINER ;  % Element containing the point
%POLYINFO.BdomRED_interp = BdomRED_interp ;
