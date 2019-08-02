problemName  = 'DenoisingEinstein';
sLoader      = SettingsLoader(problemName);
settings     = sLoader.settings;

imageProblem = DenoisingProblem(settings);
imageProblem.solve();
