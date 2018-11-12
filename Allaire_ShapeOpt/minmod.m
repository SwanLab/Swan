
function mout = minmod(phi1,phi2)
% MINMOD FUNCTION
% This function is essential for our second order
% scheme used to solve the level set equation 
% (as well as for for reinitialization).
  sphi1=sign(phi1) ;
  sphi2=sign(phi2) ;
  mout = max(0,sphi1.*sphi2).*sphi1.*min(abs(phi1),abs(phi2)) ;
end
