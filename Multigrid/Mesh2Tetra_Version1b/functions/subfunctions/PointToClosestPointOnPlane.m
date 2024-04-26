function C=PointToClosestPointOnPlane(A,B,C,P)
% C = PointToClosestPointOnPlane(A,B,C,P)
%
% inputs,
%   A,B,C : Points (x,y,z) describing a plane
%   P : A certain (x,y,z) point  
% 
% outputs,
%   C : Point on the plane which is closest to point P
%

% Calculate normal of triangle
N = cross(A-C,B-C); N=N./sqrt(sum(N.^2));
Pline = P;
Vline = N;

% Normalized plane intersection position on line [0..1]
t = dot(N, C - Pline) / (dot(N, Vline)+eps);
    
% 3D xyz Intersection Point
C =  Pline + t * Vline;
