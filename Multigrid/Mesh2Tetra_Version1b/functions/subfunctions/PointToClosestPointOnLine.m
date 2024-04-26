function [C,u]=PointToClosestPointOnLine(L1,L2,P)
% C,u = PointToClosestPointOnLine(L1,L2,P)
%
% inputs,
%   L1,L2 : Points (x,y) describing a line
%   P : A certain (x,y) point  
% 
% outputs,
%   C : Point on the line which is closest to point P
%   u : Percentage [0..1] of position on line from L1 to L2
%
u=((P(1)-L1(1))*(L2(1)-L1(1))+(P(2)-L1(2))*(L2(2)-L1(2)))/sum((L2-L1).^2);
C=L1+u*(L2-L1);
