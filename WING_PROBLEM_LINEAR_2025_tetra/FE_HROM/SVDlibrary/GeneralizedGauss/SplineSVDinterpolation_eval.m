function [f,df] = SplineSVDinterpolation_eval(xnew,Vxpoly,Uypoly,Vxpoly_der,Uypoly_der,S,DATA_SVD_INTERPOLATION)
%  This function takes as inputs the outputs of SplineSVDinterpolation_poly, 
% and a new set of points xnew, and it returns the interpolated values at
% xnew, as well as their derivatives 
% JAHO, 3-April-2020, 21th -Quarantine, COVID-19 
% -----------------------------------------------------------------------------------
npoints = size(xnew,1) ; 
f = zeros(npoints,1); 
dfx = f; 
dfy  = f ; 
for i = 1:length(S)    
  px =  Vxpoly{i} ; 
  py = Uypoly{i} ; 
  dx = Vxpoly_der{i} ; 
  dy =Uypoly_der{i}; 
  
  switch DATA_SVD_INTERPOLATION.METHOD
      case 'SPLINE'
          fx = ppval(px,xnew(:,1)) ;
          fy = ppval(py,xnew(:,2)) ;
          
          d_fx = ppval(dx,xnew(:,1)) ;
          d_fy = ppval(dy,xnew(:,2)) ;
      case 'POLYNOMIAL'
           fx = polyval(px,xnew(:,1)) ;
          fy = polyval(py,xnew(:,2)) ;
          
          d_fx = polyval(dx,xnew(:,1)) ;
          d_fy = polyval(dy,xnew(:,2)) ;
  end
  
  f = f + fx.*fy*S(i) ; 
   dfx = dfx + d_fx.*fy*S(i) ; 
  dfy = dfy + fx.*d_fy*S(i) ; 
end

df{1} = dfx; 
df{2} = dfy; 