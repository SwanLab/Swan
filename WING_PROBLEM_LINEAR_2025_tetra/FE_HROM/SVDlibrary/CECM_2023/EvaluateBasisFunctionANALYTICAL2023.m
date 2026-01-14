function [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionANALYTICAL2023(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)
%  
%   Analytical evaluation of basis functions and their gradients. 
% It also exploits the existence of an underlying FE discretization (to check whether the point is inside within the domain)
% See EvaluateBasisFunctionALL.m  


if nargin == 0
    load('tmp1_48.mat')
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
%nelem = size(VAR_SMOOTH_FE.CN,1)  ;  % Total number of FE elements

%if  isempty(POLYINFO.setElements)
    % There is no information of where the points were located in the
    % previous iteration
    [POLYINFO,ELEMENTS_CONTAINING_xNEW,PHIk_y] = LocationPointsInElement_noCN(xNEW,VAR_SMOOTH_FE,POLYINFO,PHIk_y) ;
%else
    % We have information of the location of the points in the previous
    % iteration 
 %   [POLYINFO,ELEMENTS_CONTAINING_xNEW,PHIk_y] = LocationPointsInElement_withCN(xNEW,VAR_SMOOTH_FE,POLYINFO,PHIk_y) ;
%end

if DATA.OnlyCheckIfIsInside == 0        
    [Phi_loc,Grad_loc] = LocEvalBasFunAnalytical2023(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ;    
    PHIk_y   =   Phi_loc ;
    for idim = 1:length(Grad_loc)        
        dPHIk_y{idim} = Grad_loc{idim};
    end    
end

POLYINFO.setElements = ELEMENTS_CONTAINING_xNEW ; %(inew ) = elemCONTAINER ;  % Element containing the point
 
