function main
close all
problemNames  = {'DenoisingEinstein';'AcceleratedDenoisingEinstein'};

for iproblem = 1:length(problemNames)
    pName = problemNames{iproblem};
    sLoader      = SettingsLoader(pName);
    settings     = sLoader.settings;
    denoisingProblem = DenoisingProblem(settings);
    denoisingProblem.solve();
    pD = denoisingProblem.plottingData;
    figure(1)
    hold on
    plot(pD.cost)
    figure(2)
    hold on
    plot(pD.dualGap)
end

end
