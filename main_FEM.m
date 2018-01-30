clc
clear variables
a=addpath(genpath('./FEM/'));
b=addpath(genpath('./Input/'));


%% Steps
% 1 - Run 'Main.m'
% 2 - Create object --> obj = Physical_Problem(filename);
% 3 - Preprocess    --> obj.preProcess;
% 4 - Compute       --> obj.computeVariables;
% 5 - Postprocess   --> obj.postProcess;
%% test
% run('test.m')
clear variables
%% Main.m
triangle_linear = Physical_Problem('CantileverBeam_Triangle_Linear');
triangle_linear.preProcess;
triangle_linear.computeVariables;
triangle_linear.postProcess;

fprintf('Ok\n');

% post = Postprocess_PhysicalProblem;
% gidPath = 'C:\Program Files\GiD\GiD 13.0.3\'; %write your GiD path
% files_name = triangle_linear.problemID;
% files_folder = fullfile(pwd,'Output');
% iterations = 1:1;

% output_video_name = fullfile(pwd,'StressVideo');
% post.Print_make_video_stress(gidPath,files_name,files_folder,iterations,output_video_name)


Micro_Square_Triangle = Physical_Problem_Micro('RVE_Square_Triangle');
Micro_Square_Triangle.preProcess;
Micro_Square_Triangle.computeVariables([1 0 0]);
Micro_Square_Triangle.postProcess;
Micro_Square_Triangle.computeChomog;

fprintf('Ok\n');


 
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

