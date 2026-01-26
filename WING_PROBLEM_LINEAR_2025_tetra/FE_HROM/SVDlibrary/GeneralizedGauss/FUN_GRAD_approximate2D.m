function [g,Grad_g] = FUN_GRAD_approximate2D(xNEW,DATAFITTING)
% See FUN_GRAD_approximate2D_aux.mlx

if nargin == 0
    load('tmp1.mat')
end
% ---------------------------------------------------------------------------
SingVAL = DATAFITTING.SPLINE_COEFFS_SingularValues ; 
npoints = size(xNEW,1) ; % number of points at which to evaluate the function
nfun = length(SingVAL) ;  % Number of functions to evaluate 
g = zeros(nfun,npoints) ;% Initialization   Functions to be  evaluated at points xNEW 
Grad_g = {g,g} ;  % Initialization derivatives 

  

for  iFUN = 1:nfun
    Ux = DATAFITTING.SPLINE_COEFFS_G_dim{1}{iFUN} ;  
    dUx = DATAFITTING.DER_SPLINE_COEFFS_G_dim{1}{iFUN} ;  
    Vy = DATAFITTING.SPLINE_COEFFS_G_dim{2}{iFUN} ; 
    dVy = DATAFITTING.DER_SPLINE_COEFFS_G_dim{2}{iFUN} ; 
    S = SingVAL{iFUN} ; 
    
     
    
    for  isingVAL = 1:length(S)
        ux = ppval(Ux{isingVAL},xNEW(:,1)) ;
        dux = ppval(dUx{isingVAL},xNEW(:,1)) ;
        vy = ppval(Vy{isingVAL},xNEW(:,2)) ;
        dvy = ppval(dVy{isingVAL},xNEW(:,2)) ;
        g(iFUN,:) = g(iFUN,:) + (S(isingVAL).*ux.*vy)' ;
        idim = 1; 
        Grad_g{idim}(iFUN,:) = Grad_g{idim}(iFUN,:) + (S(isingVAL).*dux.*vy)' ;
         idim = 2; 
        Grad_g{idim}(iFUN,:) = Grad_g{idim}(iFUN,:) + (S(isingVAL).*ux.*dvy)' ;        
    end
    
      
    
end

 