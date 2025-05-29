% MLX files
run('Tutorial00_Mesh.mlx');
run('Tutorial02_FEMElasticity.mlx');
run('Tutorial03_LevelSet.mlx');
run('Tutorial03p2_UnfittedMeshFunction.mlx');
run('Tutorial04_Filters.mlx');
run('Tutorial08_Projectors.mlx');
run('Tutorial09_ProjectorsQuadrilaterals.mlx');
run('Tutorial10_AcademicProblem.mlx');

% M files
TopOptTestTutorial();
TopOptTestTutorial3DDensity();
TopOptTestTutorialBoundFormulation();
TopOptTestTutorialDensityNullSpace();
TopOptTestTutorialGiD();
TopOptTestTutorialGlobalLengthScaleControl();
TopOptTestTutorialLevelSetNullSpace();
TopOptTestTutorialLSPerimeter();
TopOptTestTutorialMicro();
TopOptTestTutorialWithGiD();
%TopOptViaHomogenizationTutorial(); % ALEX
Tutorial02FEMElasticity();
Tutorial02p2FEMElasticityMicro();
Tutorial11Monitoring();
TutorialRemeshing();
ChomogNetworkTutorial();
TutorialXXPhaseFieldCase;
TutorialXXPhaseFieldHomogenization;
close all;

% Output:
clear;
close all;
clc;
fprintf('Congratulations, all tutorials are working!');
fprintf('\n');