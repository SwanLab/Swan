% InitTraining Data
clc; clear;

%% INPUTS

% MESH
s.meshSource = 'Generated';     % Generated / Imported
s.meshFile   = 'DEF_Q4auxL_1';  % Just for imported mesh
s.ndim       = 2;
s.nelem      = [40 40];         % Just for generated mesh          
s.bounds     = [3, 2];          % Just for generated mesh

% GEOMETRY - Just for generated mesh
gPar.type           = 'Auxetic';
gPar.length         = 2;
gPar.height         = 2;
gPar.xCoorCenter    = 1.5;
gPar.yCoorCenter    = 1;
gPar.theta          = 58;
gPar.beta           = 64;
gPar.thickness      = 0.3; 
s.geometry          = GeometricalFunction(gPar);

% MATERIAL
s.inclusionType = 'Void';     % Void / Weak
s.young         = 1;          
s.poisson       = 1/3;        

% TRAINING
s.method        = 'Multiscale';     % Multiscale / EIFEM
s.nSubdomains   = [5 5];
s.domainIndices = [3 3];           
s.basisType     = 'Discontinuous';  % Continuous / Discontinuous
s.coarseOrder   = 1;
s.dirac         = true;

%% COMPUTE
training = TrainingPreconditioning(s);
results  = training.compute();
