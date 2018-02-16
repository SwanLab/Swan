phi = pi/4;
theta = 0;
%Shear

% lambda = 13.61:0.05:20;
% tipo = 'MIN_INV_STIFF';
% for ilambda = 1:length(lambda) 
% MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL',tipo);
% end
% 
% lambda = 0.01:0.05:20;
% tipo = 'MIN_MINUS_STIFF';
% for ilambda = 1:length(lambda) 
% MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL',tipo);
% end


lambda = 2.46:0.05:20;
for ilambda = 1:length(lambda)   
MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL','MIN_STIFF_INV');
end


%Horizontal
phi = pi/2;
theta = pi/2;

lambda = 0.01:0.05:20;
for ilambda = 1:length(lambda) 
MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL','MIN_MINUS_STIFF');
MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL','MIN_STIFF_INV');
MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL','MIN_INV_STIFF');
end



%Bulk
phi = pi/4;
theta = pi/2;

lambda = 0.01:0.05:30;
for ilambda = 1:length(lambda) 
MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL','MIN_MINUS_STIFF');
MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL','MIN_STIFF_INV');
MAIN_OPT3(0.5,lambda(ilambda),0,phi,theta,'LINEAL','MIN_INV_STIFF');
end