% Swan blocks

% Settings reader
fileName               = 'test_arturo';
s                      = Settings(fileName);
s.warningHoleBC        = false;
s.printIncrementalIter = false;
s.printChangingFilter  = false;
s.printing             = false;
translator             = SettingsTranslator();
translator.translate(s);
fileName               = translator.fileName;
settings               = SettingsTopOptProblem(fileName);

% MESH (gid approach)
file       = 'Arturo';
a.fileName = file;
data       = FemDataContainer(a);
Mesh       = data.mesh;

% Design variable
d                               = settings.designVarSettings;
d.mesh                          = Mesh;
d.scalarProductSettings.epsilon = Mesh.computeMeanCellSize;
d.scalarProductSettings.mesh    = Mesh;
rho                             = DesignVariable.create(d);
rhoElem                         = squeeze(rho.computeVolumeFraction());

% Material
homogSettings            = settings.homogenizedVarComputerSettings;
m.constitutiveProperties = homogSettings.constitutiveProperties;
m.typeOfMaterial         = homogSettings.typeOfMaterial;
m.interpolation          = homogSettings.interpolation;
m.nElem                  = homogSettings.nelem;
m.dim                    = homogSettings.dim;
Cp                       = MaterialInterpolation.create(m);
mp                       = Cp.computeMatProp(rhoElem);
Mat                      = data.material;
Mat.compute(mp);

% FEM
f               = settings.problemData.femData;
physicalProblem = FEM.create(f);

% Shape functionals
homogVar                 = HomogenizedVarComputer.create(homogSettings);
c.homogVarComputer       = homogVar;
c.targetParameters       = [];
c.targetParameters.Vfrac = s.Vfrac_final;
c.femSettings.fileName   = file;
c.designVariable         = rho;
c.type                   = settings.constraintSettings.shapeFuncSettings{1,1}.type;
c.filterParams           = settings.constraintSettings.shapeFuncSettings{1,1}.filterParams;
compliance               = ShapeFunctional.create(c);

c.type         = settings.constraintSettings.shapeFuncSettings{2,1}.type;
c.filterParams = settings.constraintSettings.shapeFuncSettings{2,1}.filterParams;
volume         = ShapeFunctional.create(c);

% Computing function value and gradient
compliance.computeFunctionAndGradient();
volume.computeFunctionAndGradient();