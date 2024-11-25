function Lecture2
t0   = 0;
tmax = 1;
nT = 100;
h = (tmax - t0)/nT;
t = linspace(0,1,nT);

y0 = 1;
y(:,1) = y0;
dy = @(t,y) computeDerivative(t,y);
for it = 1:nT-1
    y(:,it+1) = forwardEuler(t(it),y(:,it),h,dy);   
end
plot(t,y)
end

function dy = computeDerivative(t,y)
dy = y*y*t;
end


function yt1 = forwardEuler(t,yt,h,dy)
f = dy(t,yt);
yt1 = yt + h*f;
end