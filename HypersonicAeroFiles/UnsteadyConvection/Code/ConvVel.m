function [Conv] = ConvVel(X)
% Conv = ConvVel(X)
% Given a matrix of nodal coordinates X, 
% this functions provides a nodal velocity field Conv
%



   %----- v(x,y) = (-y,x)
   Conv = zeros(size(X));
   Conv(:,1) =-X(:,2);
   Conv(:,2) = X(:,1);


