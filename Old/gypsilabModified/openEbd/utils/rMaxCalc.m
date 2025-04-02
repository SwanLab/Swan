function[rM] = rMaxCalc(X,Y)
dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
rX   = dist(X,mean(X,1));
rY   = dist(Y,mean(Y,1));
rXY  = dist(mean(X,1),mean(Y,1));
rM = rXY + rX  + rY;
end