function ENIAfirstLectureLearnMatlab()
nSteps = 10;
t = createTimeStepping(nSteps);
theta = t.*t;
plot(t,theta,'-+');
end

function t = createTimeStepping(nSteps)
t = linspace(0,1,nSteps);
end