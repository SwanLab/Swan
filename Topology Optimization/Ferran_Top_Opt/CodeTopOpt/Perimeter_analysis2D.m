function  Perimeter_analysis2D

g = [2 2 2 2;1 1 0 0;1 0 0 1; 0 1 0 0; 1 1 1 0;1 1 0 1;0 0 1 0; 0 0 0 0;0 0 0 0;0 0 0 0];


n_eps = 6;
epsilon = 0.5.^[0:n_eps]';

he = ones(n_eps+1,1)*6.25*1e-3;
he2 = [2.5*1e-2  2.5*1e-2  2.5*1e-2 1.25*1e-2 6.25*1e-3 3.125*1e-3 1.5625*1e-3 7.8125e-04 3.9063e-04 1.9531e-04]';

%he = 0.5.^[8:n_eps+8];

%figure(3)
%plot(epsilon,'-+')
%hold on
%plot(1:length(epsilon),h*ones(length(epsilon),1),'r')
mesh = [];


etas(:,1) = ones(n_eps+1,1);
etas(:,2) = 5*ones(n_eps+1,1);
etas(:,3) = 10*ones(n_eps+1,1);
etas(:,4) = epsilon./he2(1:n_eps+1);
etas(:,5) = epsilon./he;

remesh = ones(size(etas));
remesh(1,:) = zeros(1,size(etas,2));
for i = 1:size(etas,2)
hs(:,i) = epsilon./etas(:,i);

for j = 1:size(etas,1)-1
    if hs(j+1,i) == hs(j,i)
        remesh(j+1,i) = 0;
    end
end


%remesh = 
end




for ieta = 1:size(etas,2)
    
    [p,e,t] = initmesh(g(:,1:4),'hmax',hs(1,ieta), 'jiggle','mean');
    mesh.g = g; mesh.p = p; mesh.e = e; mesh.t = t;
for iepsilon = 1:n_eps+1
    
    [K,M,Nodes,mesh] = mass_and_Stiff(mesh,remesh(iepsilon,ieta));
    N(iepsilon,ieta) = size(Nodes,2);
    txi = compute_caracteristic(Nodes);
    
    eps = epsilon(iepsilon);
    % Compute ve
    ve = (eps*eps*K + M)\(M*txi);
    Pe(iepsilon) = (ve)'*M*(ve - 2*txi);
    vol = compute_volum(txi,M);
    Per(iepsilon,ieta) = 4/eps*[Pe(iepsilon) + vol];
    figure(1)
    plot(Per,'-+')
    hold on
    plot([1:n_eps+1],ones(1,n_eps+1)*3*pi/4,'-k')
    hold off
    figure(2)
    plot3(Nodes(1,:),Nodes(2,:),txi,'+')
    hold on
    plot3(Nodes(1,:),Nodes(2,:),ve,'r+')
    hold off
    
    
    index = abs(Nodes(2,:)-0.5) < 0.1*hs(iepsilon,ieta);
    figure(4)
    plot(Nodes(1,index),[txi(index),ve(index)],'+') 

    figure(5); clf; pdeplot(mesh.p,mesh.e,mesh.t,'xydata',txi,'xystyle','flat','colormap','gray','xygrid','off','colorbar','off'); axis image; axis off; 
    %     [x,y] = meshgrid(Nodes(1,:),Nodes(2,:));
%     NodesPlot = [x(:)'; y(:)'];
%     z = compute_caracteristic(NodesPl4)ot);
%     zplot = reshape(z,size(x,1),size(x,2));
%     tri = delaunay(x,y);
%     trisurf(tri,x,y,zplot)
    
    %axis([-0.5 0.5 0 1])

end

end
%plot(epsilon,Per,'+')








end

function txi = compute_caracteristic(Nodes)
txi = zeros(size(Nodes,2),1);
center = [0.5 0.5];
txi((Nodes(1,:)-center(1)).^2 + (Nodes(2,:)-center(2)).^2 <= 0.375^2) = 1;
% center = [0.25 0.75];
% txi((Nodes(1,:)-center(1)).^2 + (Nodes(2,:)-center(2)).^2 <= 0.1^2) = 0;
% center = [0.75 0.25];
% txi((Nodes(1,:)-center(1)).^2 + (Nodes(2,:)-center(2)).^2 <= 0.1^2) = 0;
% center = [0.75 0.75];
% txi((Nodes(1,:)-center(1)).^2 + (Nodes(2,:)-center(2)).^2 <= 0.1^2) = 0;

end

function vol = compute_volum(txi,M)
vol = txi'*M*txi; 
end




function [K,M,Nodes,mesh] = mass_and_Stiff(mesh,remesh)
g = [2 2 2 2;1 1 0 0;1 0 0 1; 0 1 0 0; 1 1 1 0;1 1 0 1;0 0 1 0; 0 0 0 0;0 0 0 0;0 0 0 0];
% polig=[ 2 6  0.5  0.25 -0.25 -0.5 -0.25 0.25 0 yb yb 0 -yb -yb]';
% gd = [polig];
% ns = char('polig')';
% sf = '(polig)';
% [g,bt] = decsg(gd,sf,ns); % combines the basic shapes
p = mesh.p; g = mesh.g; e = mesh.e; t = mesh.t;
if remesh == 1
[p,e,t] = refinemesh(g,p,e,t,'regular');
end
%[p,e,t] = refinemesh(g,p,e,t,'regular');
%[p,e,t] = refinemesh(g,p,e,t,'regular');

% figure(1); clf; pdemesh(p,e,t); axis image; axis off;
mesh.g = g; mesh.p = p; mesh.e = e; mesh.t = t;
c = 1;a = 1;f = 1;
[K,M,~] = assema(p,t,c,a,f);
Nodes = p;

% applyBoundaryCondition(model,'Edge',1:4,'u',0);
% generateMesh(model,'Hmax',h);
% [K,M,F] = assema(model,c,a,f);
% Nodes = model.Mesh.Nodes;
%[K,M,F,Q,G,H,R] = assempde(model,c,a,f);
%u = assempde(K,M,F,Q,G,H,R);
%pdeplot(model,'xydata',u,'zdata',u)
end