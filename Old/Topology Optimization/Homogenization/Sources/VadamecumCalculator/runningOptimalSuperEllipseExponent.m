function runningOptimalSuperEllipseExponent

s.samplePoints = createSamplePoints();
s.fileName = 'OptimalSuperEllipseExponentDataFromFixedRho';
exponentComputer = OptimalExponentComputer(s);
exponentComputer.compute();

end

function sample = createSamplePoints()
s.type = 'FromMxMy';
sample = SamplePointsCreatorForOptimalExponentComputer.create(s);

s.type = 'FromFixedRho';
s.rho0 = 0.8;
s.psi = pi/4;
sample = SamplePointsCreatorForOptimalExponentComputer.create(s);

end