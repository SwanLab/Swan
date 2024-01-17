function [u,A,b] = fempoi(p,t,e)
% p is the locations (x,y) of the node points
% t is the triangulation connectivity
% e is the index wrt p of the boundary nodes
%
%       Square, Dirichlet left/bottom:
%
%pv = [0,0; 1,0; 1,1; 0,1; 0,0];
%[p,t,e] = pmesh(pv, 0.15, 0);
%e = e(p(e,1) < 1e-6 | p(e,2) < 1e-6);
%u = fempoi(p,t,e);
%tplot(p,t,u)
%
%       Circle, all Dirichlet:
%
%n = 32; phi = 2*pi*(0:n)'/n;
%pv = [cos(phi), sin(phi)];
%[p,t,e] = pmesh(pv, 2*pi/n, 0);
%u = fempoi(p,t,e);
%tplot(p,t,u)

A = sparse(zeros(size(p,1),size(p,1))); b = sparse(zeros(size(p,1),1));
for i = 1:size(t,1)
    [Ak, bk] = assembly([p(t(i,:),1),p(t(i,:),2)]);
    A(t(i,:),t(i,:)) = A(t(i,:),t(i,:)) + Ak ;
    b(t(i,:)) = b(t(i,:)) + bk;
    %spy(A);
end

for j = 1:size(e,1)
    A(e(j), 1:size(A,1)) = 0;
    A(1:size(A,1), e(j)) = 0;
    b(e(j)) = 0;
    A(e(j), e(j)) = 1;  
end
    u = full(A\b);
end
