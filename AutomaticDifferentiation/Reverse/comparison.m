%% REVERSE vs FORWARD

vals = linspace(1,300,100000);
numOfVars = 1000:1000:100000;

tReverse = zeros(1,length(numOfVars));
tForward = tReverse;
for n = 1:length(numOfVars)
    [t, gradR] = computeReverseVarEvaluation(vals(1:numOfVars(n)));
    tReverse(n) = t;
    [t, gradF] = computeForwardVarEvaluation(vals(1:numOfVars(n)));
    tForward(n) = t;
end
figure()
plot(numOfVars,tForward,numOfVars,tReverse)
xlabel('Number of variables');
ylabel('Time (s)')
legend('Forward','Reverse')

