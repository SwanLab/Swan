function ENIAfirstLectureLearnMatlab()
nSteps = 10;
t = createTimeStepping(nSteps);
theta = t.*t;
plot(t,theta,'-+');
dtheta = t.*t.*t;
hold on
plot(t,dtheta,'-+');
ddtheta = log(t).*t;
figure
plot(t,ddtheta,'-+');
end

function t = createTimeStepping(nSteps)
t = linspace(0,1,nSteps);
end