function dispMACRO = CoarseDisplace(strainINP,COORref,DATA,h,dispREF)




ex =  strainINP(1) ;
ey = strainINP(2) ;
exy = 0.5*strainINP(3) ;
exz = 0.5*strainINP(8) ;
eyz = 0.5*strainINP(7) ;
kx =  strainINP(4) ;
ky = strainINP(5) ;
kxy = 0.5*strainINP(6) ;

%
MACRODEF_0 = [ex exy exz
    exy ey eyz ] ;
MACRODEF_1 = [kx kxy
    kxy ky ] ;

% u,v
z = COORref(3,:) ;
zCOORref = bsxfun(@times,COORref(1:2,:)',z') ;
dispMACROuv =  MACRODEF_0*COORref(1:3,:)   + MACRODEF_1*zCOORref';



%%%%% 3 order theory
if   DATA.ORDER_SHEAR_THEORY == 3 % | DATA.SPECIAL_BC_SHEAR_ON == 13;
    MACRODEF_0 = [ex exy 1.5*exz
        exy ey 1.5*eyz ] ;
    MACRODEF_1 = [kx kxy
        kxy ky ] ;
    % u,v
    z = COORref(3,:) ;
    zCOORref = bsxfun(@times,COORref(1:2,:)',z') ;
    dispMACROuv =  MACRODEF_0*COORref(1:3,:)   + MACRODEF_1*zCOORref';
    dispWARPING_xz = -4*z.^3/3/h^2*[exz ]*(3/2) ;
    dispWARPING_yz = -4*z.^3/3/h^2*[ eyz]*(3/2) ;
    dispMACROuv = dispMACROuv + 2*[dispWARPING_xz; dispWARPING_yz ];
    wL =  1.5*exz*COORref(1,:)+1.5*eyz*COORref(2,:) - kx*COORref(1,:).^2 ...
        - ky*COORref(2,:).^2 -2*kxy*COORref(1,:).*COORref(2,:);
else
    MACRODEF_0 = [ex exy exz
        exy ey eyz ] ;
    MACRODEF_1 = [kx kxy
        kxy ky ] ;
    % u,v
    z = COORref(3,:) ;
    zCOORref = bsxfun(@times,COORref(1:2,:)',z') ;
    dispMACROuv =  MACRODEF_0*COORref(1:3,:)   + MACRODEF_1*zCOORref';
    wL =  exz*COORref(1,:)+eyz*COORref(2,:) - kx*COORref(1,:).^2 ...
        - ky*COORref(2,:).^2 -2*kxy*COORref(1,:).*COORref(2,:);
end


%%%%%
% w



dispMACRO = [dispMACROuv;wL] ;

if ~isempty(dispREF)
    dispMACRO =   bsxfun(@minus,dispMACRO,dispREF) ;
end