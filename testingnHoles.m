% Testing nHoles

clear
clc

fileName = 'test_micro3d';
s = Settings(fileName);
s.warningHoleBC = false;
s.printIncrementalIter = false;
s.printChangingFilter = false;
s.printing = false;
translator = SettingsTranslator();
translator.translate(s);
fileName = translator.fileName;
settings  = SettingsTopOptProblem(fileName);
settings.designVarSettings.creatorSettings.type = 'holes';
m = settings.problemData.femData.mesh;

l.mesh = m;
l.type = 'LevelSet';
l.initialCase = 'holes';
l.creatorSettings = settings.designVarSettings.creatorSettings;
LS = LevelSet(l);
p1.fValues = LS.value;
p1.mesh    = m;
Fun        = P1Function(p1);

p1.filename = 'examplenHoles';
p1.type = 'Paraview';
Fun.print(p1);