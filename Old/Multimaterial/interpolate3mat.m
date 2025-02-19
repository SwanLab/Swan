function interpolate3mat

X(:,1) = [1 0 0];
X(:,2) = [0 1 0];
X(:,3) = [0 0 1];

xi1Init = 0.4;
xi2Init = 0.01;

alpha1Init = 17;
alpha2Init = 2;
alpha3Init = 4;

rho(1,1) = xi1Init;
rho(2,1) = xi2Init;
rho(3,1) = 1 - rho(1) - rho(2);

gamma(1,1) = alpha1Init;
gamma(2,1) = alpha2Init;
gamma(3,1) = alpha3Init;

figure('color','w')
hold on

for k=1:15
    
       plotTriangle(X)

    eta12 = rho(1) / (rho(1)+rho(2));
    eta21 = rho(2) / (rho(1)+rho(2));
    
    eta13 = rho(1) / (rho(1)+rho(3));
    eta31 = rho(3) / (rho(1)+rho(3));
    
    eta23 = rho(2) / (rho(2)+rho(3));
    eta32 = rho(3) / (rho(2)+rho(3));
    
    B = [0 eta23 eta32;...
         eta13 0 eta31; ...
         eta12 eta21 0];
     
    %interpolation
    gamma = B*gamma;
    
    
    X = B*X
    

    
    H = 0.5*[0 1 1;...
             1 0 1;...
             1 1 0];
    
    rho = H*rho
    

end
% 
% alpha1
% alpha2
% alpha3

gamma

%linear interpolation should give the same as the following:
xi1Init * alpha1Init + xi2Init * alpha2Init + (1-xi1Init-xi2Init) * alpha3Init

end

function plotTriangle(p)

p1=p(1,:);
p2=p(2,:);
p3=p(3,:);
% Plot trianle using 'patch'

h=patch('Faces',1:3,'Vertices',[p1;p2;p3]);
set(h,'FaceColor','r','EdgeColor','k','LineWidth',2,'FaceAlpha',0.5)
axis equal vis3d
view([45 25 45])
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)
zlabel('z','FontSize',20)


end

