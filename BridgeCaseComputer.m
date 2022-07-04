function BridgeCaseComputer 
    fileName = 'test_bridg_multipleLoads';
    settingsTopOpt = SettingsTopOptProblem(fileName);            
    m = settingsTopOpt.problemData.femData.mesh;
    x = m.coord(:,1);
    y = m.coord(:,2);
    isNodeFixedUp1  = y >= 1.9;
    isNodeFixedUp2  = y >= 0.1 & (x <= 1 | x>= 9);
    isNodeFixed     = isNodeFixedUp1 | isNodeFixedUp2;
    settingsTopOpt.designVarSettings.isFixed.values = 1*ones(sum(isNodeFixed),1);
    settingsTopOpt.designVarSettings.isFixed.nodes = find(isNodeFixed);
    topOptProblem = TopOpt_Problem(settingsTopOpt);
    topOptProblem.computeVariables;
%     topOptProblem.postProcess;
end
