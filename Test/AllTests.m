FemTestsSuite;
MicroShapeTestsSuite; %OK
DiffReactTestsSuite; %OK
UnfittedIntegrationTestsSuite; %OK
VectorizedTriangulationTestsSuite; %OK
close all;
TopOptTestsSuite;
MultimaterialTestsSuite;
ReadingFilesTestsSuite;
%PlottingTestsSuite; % Funciona quan es generen els .mat
close all;
%HomogenizationTestsSuite; % No funciona (issue Swan)
ImageProcessingTestsSuite;
%PostProcTestsSuite; S'ha de fer des de 0 i pensar en com (ara ignore it)
ProjectorsTestsSuite;
RemeshingTestsSuite;
%DehomogenizationTestsSuite; % No funciona (issue Swan)
close all;
AcademicTestsSuite; %OK
close all;
BoundaryCondTestsSuite;
GeomFunTestsSuite; %OK
close all;
DomainFunTestsSuite; %OK
AlgebraicOperationsTestsSuite; %OK
PhaseFieldTestsSuite;
ContinuumDamageTestsSuite;
HyperelasticityTestsSuite;
close all;
