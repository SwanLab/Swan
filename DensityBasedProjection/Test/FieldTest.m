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
load('C:\Users\artur\Documents\GitHub\SWAM\Swan\DensityBasedProjection\Test\Data\ResultsData3Iterations.mat')
elseif s.iterations == 5
load('C:\Users\artur\Documents\GitHub\SWAM\Swan\DensityBasedProjection\Test\Data\ResultsData5Iterations.mat')
else 
warning('No test Data for the current iterations')
end
%% Create the objects
Test1 = ComplianceRobustComputer(s);
Test1.compute();

%% Validator
if results.projectedField.E == Test1.projectedField.E
    %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
    disp('Projected Field E_|OK!|')
else
    warning('Error in Projected Field E')
end

if results.projectedField.I == Test1.projectedField.I
    %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
    disp('Projected Field I_|OK!|')
else
    warning('Error in Projected Field I')
end
if results.projectedField.D == Test1.projectedField.D
    %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
    disp('Projected Field D_|OK!|')
else
    warning('Error in Projected Field D')
end



