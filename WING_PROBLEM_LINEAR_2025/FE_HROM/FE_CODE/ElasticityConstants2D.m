function [Ex,Ey,nuXY,nuYX,gXY] = ElasticityConstants2D(Celas)



D = inv(Celas) ;

Ex = 1/D(1,1) ;
Ey = 1/D(2,2) ;


 

nuXY = -D(1,2)*Ey ;
 

nuYX =  -D(2,1)*Ex ;
 
 
gXY = 1/D(3,3) ;

disp('----------------------')
disp('----------------------')
disp('Young moduli')
disp('-----------------------')
disp(['E_x = ',num2str(Ex)]) ;
disp(['E_y = ',num2str(Ey)]) ;
disp('Poisson ratios')
disp('-----------------------')
disp(['nu_xy = ',num2str(nuXY)]) ;
disp(['nu_yx = ',num2str(nuYX)]) ;
disp('Shear moduli')
disp('-----------------------')
disp(['G_xy = ',num2str(gXY)]) ;