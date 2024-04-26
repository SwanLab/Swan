function [inter,inter_x,inter_y,inter_z]=LineTriangleIntersection(A,B,C,P1,P2,ignore_corners)
if(nargin<6), ignore_corners=false; end

% http://www.cs.brown.edu/~scd/facts.html

% Calculate normal of triangle
N = cross(A-C,B-C); N=N./sqrt(sum(N.^2));
Pline = P1; Vline = P2-P1;
% Normalized plane intersection position on line [0..1]
t = dot(N, C - Pline) / (dot(N, Vline)+eps);

% 3D xyz Intersection Point
P =  Pline + t * Vline;

inter_x=P(1); inter_y=P(2); inter_z=P(3);

% If no intersection between points defining the line, return
if((t<0)||(t>1)), inter=false; return; end

% Drop x,y or z to get a 2D triangle intersection problem
[temp,i]=max(N);
P_2D=P; P_2D(i)=[]; 
A_2D=A; A_2D(i)=[]; 
B_2D=B; B_2D(i)=[];
C_2D=C; C_2D(i)=[];

% Boundary box
mi=[min([A_2D(1) B_2D(1) C_2D(1)]) min([A_2D(2) B_2D(2) C_2D(2)])];
ma=[max([A_2D(1) B_2D(1) C_2D(1)]) max([A_2D(2) B_2D(2) C_2D(2)])];

% Check if the insection point is inside the face
inter=CheckInsideFace(A_2D,B_2D,C_2D,P_2D(1),P_2D(2),mi,ma,ignore_corners);
