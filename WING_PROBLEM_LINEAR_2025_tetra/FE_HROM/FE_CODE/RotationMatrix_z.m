function T = RotationMatrix_z(R)





cos_a = R(1,1) ;
sin_a = R(2,1) ;

cosTc = (cos_a)^2 ;
sinTc = (sin_a)^2 ;
sinTd = 2*cos_a*sin_a ; %  (sin(2*a)) ;
cosTd = cosTc -sinTc ; % cos(2*a) ;


cosT = cos_a ;
sinT = sin_a ;

T = [cosTc  sinTc    0      0      0      -sinTd
    sinTc  cosTc    0      0      0      sinTd
    0      0       1      0      0        0
    0      0     0      cosT      sinT           0
    0      0        0       -sinT      cosT          0
    0.5*(sinTd)     -0.5*(sinTd)        0       0    0          cosTd         ] ;


if size(R,2) == 2
    % See /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/...
    % APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/FE_CODE/QtransfBvect_2Dps.m
    INDICES = [1 2 6 3] ; 
    T = T(INDICES,INDICES) ; 
end