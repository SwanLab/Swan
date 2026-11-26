% InitTraining Data
clc; clear;

%% INPUTS

% MESH
s.meshSource = 'Generated';   % 'Generated' / 'Imported'
s.meshFile   = 'FileName';    % Just for imported mesh
s.ndim       = 3;
s.nelem      = [10 10 10];    % Just for generated mesh  
s.bounds     = [2 2 2];       % Just for generated mesh

% GEOMETRY - Just for generated mesh
gPar.type         = 'SphereInclusion';
gPar.radius       = 0.5;
gPar.xCoorCenter  = 1;
gPar.yCoorCenter  = 1;
gPar.zCoorCenter  = 1;
s.geometry        = GeometricalFunction(gPar);


% MATERIAL
s.inclusionType = 'Void';     % Void / Weak
s.young         = 1;          
s.poisson       = 1/3;        

% TRAINING
s.method        = 'Multiscale';         % Multiscale / EIFEM
s.nSubdomains   = [5 5 5];
s.domainIndices = [3 3 3];           
s.basisType     = 'Discontinuous';      % Continuous / Discontinuous
s.coarseOrder   = 1;
s.dirac         = true;

%% COMPUTE
training = TrainingPreconditioning(s);
results  = training.compute();
