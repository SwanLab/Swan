function [g,Grad_g] = FUN_GRAD_approximate3D(xNEW,DATAFITTING)
% See FUN_GRAD_approximate3D_aux.mlx

if nargin == 0
    load('tmp1.mat')
end

% DATAFITTING = 
% 
%   struct with fields:
% 
% %                SPLINE_COEFFS_G_X: {64×1 cell}
%                   SPLINE_COEFFS_G_Y: {64×1 cell}
%                   SPLINE_COEFFS_G_Z: {64×1 cell}
%               DER_SPLINE_COEFFS_G_X: {64×1 cell}
%               DER_SPLINE_COEFFS_G_Y: {64×1 cell}
%               DER_SPLINE_COEFFS_G_Z: {64×1 cell}
%      SPLINE_COEFFS_SingularValues_Z: {64×1 cell}
%     SPLINE_COEFFS_SingularValues_XY: {64×1 cell}



% ---------------------------------------------------------------------------
 
npoints = size(xNEW,1) ; % number of points at which to evaluate the function
nfun = length(DATAFITTING.SPLINE_COEFFS_SingularValues_Z) ;  % Number of functions to evaluate 
g = zeros(nfun,npoints) ;% Initialization   Functions to be  evaluated at points xNEW 
Grad_g = {g,g,g} ;  % Initialization derivatives 

  

for  iFUN = 1:nfun
    
    spline_wz = DATAFITTING.SPLINE_COEFFS_G_Z{iFUN} ;  
    der_spline_wz = DATAFITTING.DER_SPLINE_COEFFS_G_Z{iFUN} ;       
    Sz = DATAFITTING.SPLINE_COEFFS_SingularValues_Z{iFUN} ; 
    
   
    
    for  j = 1:length(Sz)
        wz = ppval(spline_wz{j},xNEW(:,3)) ; 
        dwz = ppval(der_spline_wz{j},xNEW(:,3)) ; 
        %----------------------------------------
       % gxy = zeros(npoints,1) ; 
        S_xy = DATAFITTING.SPLINE_COEFFS_SingularValues_XY{iFUN}{j} ;
        Ux = DATAFITTING.SPLINE_COEFFS_G_X{iFUN}{j} ; 
        dUx = DATAFITTING.DER_SPLINE_COEFFS_G_X{iFUN}{j} ;
        
        Vy = DATAFITTING.SPLINE_COEFFS_G_Y{iFUN}{j} ; 
        dVy = DATAFITTING.DER_SPLINE_COEFFS_G_Y{iFUN}{j} ;
        for k=1:length(S_xy)
            ux = ppval(Ux{k},xNEW(:,1)) ;
            dux = ppval(dUx{k},xNEW(:,1)) ;
            
            vy = ppval(Vy{k},xNEW(:,2)) ;
            dvy = ppval(dVy{k},xNEW(:,2)) ;
            g(iFUN,:) = g(iFUN,:) + Sz(j)*S_xy(k)*(wz.*ux.*vy)' ;
            idim = 1 ; 
            Grad_g{idim}(iFUN,:) =Grad_g{idim}(iFUN,:)+ Sz(j)*S_xy(k)*(wz.*dux.*vy)' ;
            idim = 2 ; 
            Grad_g{idim}(iFUN,:) =Grad_g{idim}(iFUN,:)+ Sz(j)*S_xy(k)*(wz.*ux.*dvy)' ;
            idim = 3 ; 
            Grad_g{idim}(iFUN,:) =Grad_g{idim}(iFUN,:)+ Sz(j)*S_xy(k)*(dwz.*ux.*vy)' ;
            
        end
        
           
    end
    
      
    
end

 