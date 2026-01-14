function [Ex,Ey,Ez,nuXY,nuXZ,nuYX,nuZY,gYZ,gXZ,gXY] = ElasticityConstants(Celas)



D = inv(Celas) ;

Ex = 1/D(1,1) ;
Ey = 1/D(2,2) ;
Ez = 1/D(3,3) ;

nuXY = -D(1,2)*Ey ;
nuXZ = -D(1,3)*Ez ;

nuYX =  -D(2,1)*Ex ;
nuYZ = -D(2,3)*Ez ;

nuZX =  -D(3,1)*Ex ;
nuZY = -D(3,2)*Ey ;

gYZ = 1/D(4,4) ;
gXZ = 1/D(5,5) ;
gXY = 1/D(6,6) ;

disp('----------------------')
disp('----------------------')
disp('Young moduli')
disp('-----------------------')
disp(['E_x = ',num2str(Ex)]) ;
disp(['E_y = ',num2str(Ey)]) ;
disp(['E_z = ',num2str(Ez)]) ;
disp('Poisson ratios')
disp('-----------------------')
disp(['nu_xy = ',num2str(nuXY)]) ;
disp(['nu_xz = ',num2str(nuXZ)]) ;
disp(['nu_yx = ',num2str(nuYX)]) ;
disp(['nu_yz = ',num2str(nuYZ)]) ;
disp(['nu_zx = ',num2str(nuZX)]) ;
disp(['nu_zy = ',num2str(nuZY)]) ;
disp('Shear moduli')
disp('-----------------------')
disp(['G_yz = ',num2str(gYZ)]) ;
disp(['G_xz = ',num2str(gXZ)]) ;
disp(['G_xy = ',num2str(gXY)]) ;