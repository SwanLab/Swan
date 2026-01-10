function [a,b,c,d] = ConstantMatricesTransformation(s,l,L,GEOref,CROSS_SECTION,FIELDS)

YMAX = GEOref.(FIELDS{1}) ; 
YMIN = GEOref.(FIELDS{2}) ; 
FUNmax = CROSS_SECTION.(FIELDS{1}).FUN  ; 
FUNmin = CROSS_SECTION.(FIELDS{2}).FUN  ; 
dY = YMAX-YMIN  ;
ARGSmax = CROSS_SECTION.(FIELDS{1}).ARGS ;
ARGSmin = CROSS_SECTION.(FIELDS{2}).ARGS ;

% ----------------------------
% YMAX
ymax = feval(FUNmax,s,L,YMAX,dY,ARGSmax) ; 
ymax_0 = ymax(1) ; % for X = 0
ymax_F = ymax(2) ; % for X = XMAX
% YMIN
ymin = feval(FUNmin,s,L,YMIN,dY,ARGSmin) ; 
ymin_0 = ymin(1) ; % for X = 0
ymin_F = ymin(2) ; % for X = XMAX
% 
dy_0 = ymax_0-ymin_0 ;
dy_F = ymax_F - ymin_F ; 

c = dy_0/dY ; 
a = (ymax_0 -c*YMAX); 

d = 1/l*(dy_F/dY - c) ;
b = (ymax_F -a- (c + d*l)*YMAX)/l  ; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = FSINUSOIDAL(s,L,y0,dy,ARGS)
% ARGS(1) = Relative amplitude 
% ARGS(2) = Number of semi-periods
 
y = y0 + ARGS(2)*dy*sin(ARGS(1)*pi*s/L) ; 
 
end

 


function [y] = FCONSTANT(s,L,y0,dy,ARGS)
% ARGS(1) = Relative amplitude 
 
y = y0 + ARGS(1)*dy ;
y = y*ones(size(s)) ; 
 
end


function [y] = FLINEAR(s,L,y0,dy,ARGS)
% ARGS(1) = Relative amplitude at the beginning of the beam
% ARGS(2) = Relative amplitude at the end of the beam
% ARGS(3) = It is a maximum coordinate or minimum

AMPL_INI = ARGS(1) ;
AMPL_FIN = ARGS(2) ;
ISMAX = ARGS(3) ;

if  ISMAX ==1
y_CENT_old = y0 - dy/2 ; 

yINI = y_CENT_old + AMPL_INI*dy/2 ; 
yFIN = y_CENT_old + AMPL_FIN*dy/2 ; 

else 
    
    y_CENT_old =  y0 + dy/2  ; 

yINI = y_CENT_old - AMPL_INI*dy/2 ; 
yFIN = y_CENT_old - AMPL_FIN*dy/2 ; 

end
y = yINI + (yFIN-yINI)/L*s ;

 
end
