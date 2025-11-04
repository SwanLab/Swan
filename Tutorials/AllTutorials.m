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
Tutorial05_1_TopOpt2DDensityMacroMMA();
Tutorial05_2_TopOpt2DDensityMacroNullSpace();
Tutorial05_3_TopOpt2DLevelSetMacroNullSpace();
Tutorial05_4_TopOpt2DLevelSetPerimeter();
Tutorial05_5_TopOpt2DLevelSetInfillNullSpace();
Tutorial05_6_TopOpt2DLevelSetMacroGlobalLengthControl();
Tutorial05_7_TopOptDensityMacroGiD();
Tutorial05_8_TopOpt3DDensityMacroMMA();
Tutorial05_9_TopOpt2DDensityMicroNullSpace();
Tutorial05_10_TopOptDensityBoundFormulationMacro();
Tutorial05_11_TopOpt3DDensityMacroPython();
NullSpaceVerification();
%TopOptViaHomogenizationTutorial(); % ALEX
Tutorial02FEMElasticity();
Tutorial02p2FEMElasticityMicro();
Tutorial11Monitoring();
TutorialRemeshing();
verification_jacobian();
testGenDB();
OptimizeEllipseTutorial();
OptimizeSuperformTutorial();
TutorialXXPhaseFieldCase();
TutorialXXPhaseFieldHomogenization();
TutorialXXHyperelasticity();
TutorialXXContinuumDamage();
TutorialElasticityAMG();
close all;

% Output:
clear;
close all;
clc;
fprintf('Congratulations, all tutorials are working!');
fprintf('\n');