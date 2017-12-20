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
test_triangle_linear = Physical_Problem('CantileverBeam_Triangle_Linear');
test_triangle_linear.preProcess;
test_triangle_linear.computeVariables;
test_triangle_linear.postProcess;

post = Postprocess_PhysicalProblem;
gidPath = '/opt/GiDx64/13.0.2/'; %write your GiD path
files_name = test_triangle_linear.problemID;
files_folder = fullfile(pwd,'Output');
iterations = 1:1;

output_video_name = fullfile(pwd,'StressVideo');
post.Print_make_video_stress(gidPath,files_name,files_folder,iterations,output_video_name)



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
% 
% quadrilateral_serendipity = Physical_Problem('Cantileverbeam_Quadrilateral_Serendipity');
% quadrilateral_serendipity.preProcess
% quadrilateral_serendipity.computeVariables;
% quadrilateral_serendipity.postProcess;

