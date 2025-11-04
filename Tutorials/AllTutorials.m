run('Tutorial00_1_Mesh.mlx');
Tutorial00_2_Remeshing();

% Tutorial01_FEMThermal  -- PENDING

run('Tutorial02_FEMElasticity.mlx');
Tutorial02FEMElasticity();
Tutorial02p2FEMElasticityMicro();
Tutorial02p3ElasticityAMG();

run('Tutorial03_LevelSet.mlx');
run('Tutorial03p2_UnfittedMeshFunction.mlx');

run('Tutorial04_Filters.mlx');

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
%TopOptViaHomogenizationTutorial(); % ALEX
Tutorial05_11_TopOpt3DDensityMacroPython();

% Tutorial06_ShapeOptimization -- PENDING

Tutorial07_1_PhaseFieldCase();
Tutorial07_2_PhaseFieldHomogenization();

Tutorial08_ContinuumDamage();

run('Tutorial09_1_Projectors.mlx');
run('Tutorial09_2_ProjectorsQuadrilaterals.mlx');

run('Tutorial10_AcademicProblem.mlx');

Tutorial11Monitoring();

% Tutorial12_Dehomogenization -- PENDING

Tutorial13_Hyperelasticity();

NullSpaceVerification();

verification_jacobian();
testGenDB();
OptimizeEllipseTutorial();
OptimizeSuperformTutorial();

close all;

% Output:
clear;
close all;
clc;
fprintf('Congratulations, all tutorials are working!');
fprintf('\n');