clc
clear variables
addpath(genpath('.\FEM\'));
addpath(genpath('.\Input\'));

%% Steps
% 1 - Run 'Main.m'
% 2 - Create object --> obj = Physical_Problem(filename);
% 3 - Preprocess    --> obj.preProcess;
% 4 - Compute       --> obj.computeVariables;
% 5 - Postprocess   --> obj.postProcess;
%% test

run('test.m')
clear variables

%% Main.m
% triangle_linear = Physical_Problem('CantileverBeam_Triangle_Linear');
% % props.kappa = 1; props.mu = 0.4;
% % triangle_linear.setMatProps(props);
% triangle_linear.preProcess;
% triangle_linear.computeVariables;
% triangle_linear.postProcess;
% 
% 
% triangle_quadratic = Physical_Problem('CantileverBeam_Triangle_Quadratic');
% triangle_quadratic.preProcess;
% triangle_quadratic.computeVariables;
% triangle_quadratic.postProcess;
% 
% 
% quadrilateral_bilinear = Physical_Problem('Cantileverbeam_Quadrilateral_Bilinear');
% quadrilateral_bilinear.preProcess;
% quadrilateral_bilinear.computeVariables;
% quadrilateral_bilinear.postProcess;

quadrilateral_serendipity = Physical_Problem('Cantileverbeam_Quadrilateral_Serendipity');
quadrilateral_serendipity.preProcess
quadrilateral_serendipity.computeVariables;
quadrilateral_serendipity.postProcess;

