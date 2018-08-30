function gout = g(u1,u2,v1,v2)
% FLUX FUNCTION
% This is the numerical flux of our scheme for Hamilton Jacobi equations
  gout = sqrt( max(u1,0.).^2 + min(u2,0.).^2+ max(v1,0.).^2 + min(v2,0.).^2 ) ;
end
