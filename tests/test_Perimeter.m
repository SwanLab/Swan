clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

filenames='CantileverTriangleFine_Case_1_1_1';
settings=Settings(filenames);
mesh = Mesh_GiD(settings.filename);
settings.pdim = mesh.pdim;
settings.target_parameters.epsilon_perimeter = mesh.mean_cell_size;
settings.plotting = true;
settings.monitoring = false;
settings.initial_case = 'circle';


Monitoring = Monitoring.create(settings,mesh,settings.monitoring,settings.plotting);
design_variable_initializer = DesignVaribleInitializer.create(settings,mesh,settings.target_parameters.epsilon_perimeter);
x = design_variable_initializer.compute_initial_design();
Monitoring.plotX(x);
Perimeter = ShFunc_Perimeter(settings);
Perimeter.filter.preProcess();
Perimeter.filter.diffReacProb.preProcess();
Perimeter.computeCostAndGradient(x);
NumericalPerimeter = Perimeter.value ;
AnalyticalCirclePerimeter = 2*pi*design_variable_initializer.radius;

TestPassed = ((NumericalPerimeter - 0.3270) < 1e-4) & ((AnalyticalCirclePerimeter - 1.2566) < 1e-4);

if TestPassed
   cprintf('green','Succeed Test') 
else
   cprintf('red','Not succeed Test')  
end