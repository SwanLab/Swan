function [x0,y0]=plot_p2(sol,xnode)

numnp = size(xnode,2);
numel = (numnp-1)/2;

x0 =[];
y0 =[];

x = [-1:0.1:1];
N = [((x-1).*x/2)'  (1-x.^2)'  ((x+1).*x/2)'];

for i=1:numel
   isp = [2*i-1 2*i 2*i+1];
   delta = (xnode(2*i+1)-xnode(2*i-1))/2;
   x0 = [x0,(xnode(2*i)+x*delta)];
   y = N*sol(isp);
   y0 = [y0,y'];
end