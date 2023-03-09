%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   UNIT TESTING SCRIPT Global   %%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
%% Define the number of iterations
s.iterations = 3;

%% Load results
if s.iterations == 3
    file = fullfile("DensityBasedProjection",'Test','Data','ResultsData3Iterations.mat');    
    load(file)
elseif s.iterations == 5
    file = fullfile("DensityBasedProjection",'Test','Data','ResultsData5Iterations.mat');
    load(file)
else 
error('No test Data for the current iterations')
end
%% Create the objects
Test1 = ComplianceRobustComputer(s);
Test1.compute();

%% Validator
if results.projectedField.E == Test1.projectedField.E
    %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
    disp('Projected Field E |OK!|')
else
    warning('Error in Projected Field E')
end

if results.projectedField.I == Test1.projectedField.I
    %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
    disp('Projected Field I |OK!|')
else
    warning('Error in Projected Field I')
end
if results.projectedField.D == Test1.projectedField.D
    %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
    disp('Projected Field D |OK!|')
else
    warning('Error in Projected Field D')
end
close all 


