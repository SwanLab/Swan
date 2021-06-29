function [ X ] = uniformCircle(c,R,N)
%  X  = uniformCircle(c,R,N)
% inputs : c = [c1,c2] center of the circle, R radius, N number of points
% output : X of size Nx2, points uniformly distributed on the cricle
theta = 2*pi*rand(N,1);
X = [c(1)+R*cos(theta) c(2)+R*sin(theta)];

end

