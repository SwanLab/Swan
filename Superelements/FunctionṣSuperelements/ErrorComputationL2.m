function [DispError, StressError] = ErrorComputationL2(refTestValues,Results,meshDecomposed)
    subMesh = meshDecomposed.subMesh;
    subBoundMesh = meshDecomposed.subBoundMesh;

    refU = refTestValues.refU;
    refSigma = refTestValues.refSigma;
    refTau = refTestValues.refTau;

    u = Results.u{end};
    stress = Results.stress{1};

    % Displacement Error
    idxTip = subMesh(end).coord(:,1) == 1 & subMesh(end).coord(:,2) == 0.25;
    uTip   = full(u(idxTip,2));

    DispError.total = sqrt((refU-uTip).^2);
    DispError.relative = DispError.total/abs(refU);

    % Stress error
    idx = subBoundMesh(1,1).mesh.coord(:,1) == 0;
    Y = subBoundMesh(1,1).mesh.coord(idx,2);
    
    Sigma = [interp1(Y,stress.fValues(:,1),0.05) interp1(Y,stress.fValues(:,1),0.1)];
    Tau   = [interp1(Y,stress.fValues(:,2),0.05) interp1(Y,stress.fValues(:,2),0.1)];
   
    StressError.total = [sqrt((Sigma-refSigma).^2); sqrt((Tau-refTau).^2)];
    StressError.relative = [abs(StressError.total(1,:)./refSigma);  abs(StressError.total(2,:)./refTau)];
end