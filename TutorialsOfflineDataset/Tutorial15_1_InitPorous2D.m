% InitTraining Data
clc; clear;

%% INPUTS

% MESH
s.meshSource = 'Generated';   % Generated / Imported
s.meshFile   = 'FileName';    % Just for imported mesh
s.ndim       = 2;     
s.nelem      = [20 20];       % Just for generated mesh
s.bounds     = [2 2];         % Just for generated mesh

% GEOMETRY - Just for generated mesh
g.type        = 'CircleInclusion';
g.radius      = 0.7;
g.xCoorCenter = 1;
g.yCoorCenter = 1;
s.geometry    = GeometricalFunction(g);

% MATERIAL
s.inclusionType = 'Weak';  % Void / Weak
s.young         = 1;       
s.poisson       = 1/3;     

%TRAINING
s.method        = 'EIFEM';         % Multiscale / EIFEM
s.nSubdomains   = [5 5];
s.domainIndices = [3 3];           
s.basisType     = 'Continuous'; % Continuous / Discontinuous
s.coarseOrder   = 1;
s.dirac         = true;


%% COMPUTE
training = TrainingPreconditioning(s);
results  = training.compute();
