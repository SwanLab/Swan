function [ d ] = dist2(x0,Y)
% Evaluate the distance between a point X in R^2 and each point Y(1,:)

d = sqrt((x0(1)-Y(:,1)).^2+(x0(2)-Y(:,2)).^2);

end

