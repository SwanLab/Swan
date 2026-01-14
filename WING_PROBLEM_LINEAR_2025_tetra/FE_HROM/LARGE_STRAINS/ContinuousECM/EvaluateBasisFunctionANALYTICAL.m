function [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionANALYTICAL(xNEW,DATA,VAR_SMOOTH_FE)
%  
%   Analytical evaluation of basis functions and their gradients. 
% It also exploits the existence of an underlying FE discretization (to check whether the point is inside within the domain)
% See EvaluateBasisFunctionALL.m  


if nargin == 0
    load('tmp.mat')
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


[INDEXES_NEAR,IND_POLYG] = FindClosestNodes(xNEW,VAR_SMOOTH_FE) ; 
ELEMENTS_CONTAINING_xNEW = zeros((npoints),1) ; % Indices of the elements containing the points
%
%INDfin = 0 ;
inew = 1;
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
    
    

    
    inew = inew +1 ;
end

if DATA.OnlyCheckIfIsInside == 0
    
    %         % Determining the matrices for constructing shape functions
    %          NODESloc = small2large(elemCONTAINER,VAR_SMOOTH_FE.ngausE) ;
    %         [COEFFSpol,LelemSCALING,coorREF,POLYINFO] = CoeffsPolyShapeFunctions(elemCONTAINER,NODESloc,VAR_SMOOTH_FE,POLYINFO)  ;
    %
    %         % Evaluating the shape functions (and derivatives) at  point xLOC
    %         % ---------------------------------------------------------------------
    %         xLOC  = (xLOC-coorREF)./LelemSCALING   ;  % Transformed coordinate of the point under study
    %         factorDERIVATIVE = 1./LelemSCALING ;    % This for computing the derivatives
    %          % -------------------------------------------------------------------------------------------
    %         % Computatio  of the shape functions (recall that the "nodes" here of such functions are the Gauss points)
    %         [Pevaluate, PevaluateDER]= CoordinateMatrixPolynomial(xLOC,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
    %         N = Pevaluate*COEFFSpol ;  % Shape functions at the given points COORevaluate
    
    [Phi_loc,Grad_loc] = LocEvalBasFunAnalytical(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ;
    
    PHIk_y   =   Phi_loc ;
    for idim = 1:length(Grad_loc)
        
        dPHIk_y{idim} = Grad_loc{idim};
    end
    
end







POLYINFO.ELEMENTS_CONTAINING_xNEW = ELEMENTS_CONTAINING_xNEW ; %(inew ) = elemCONTAINER ;  % Element containing the point
%POLYINFO.BdomRED_interp = BdomRED_interp ;
