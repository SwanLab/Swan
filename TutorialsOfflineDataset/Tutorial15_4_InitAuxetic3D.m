% InitTraining Data  (STILL ONGOING)
clc; clear;

%% INPUTS

% MESH
s.meshSource = 'Imported';     %  Generated / Imported
s.meshFile   = 'DEF_por3D';    % Just for imported mesh
s.ndim       = 3;

% MATERIAL
s.inclusionType = 'Void';     % Void / Weak
s.young         = 1;          
s.poisson       = 1/3;        

% TRAINING
s.method        = 'EIFEM';     % Multiscale / EIFEM
s.nSubdomains   = [3 3 1];
s.domainIndices = [2 2 1];           
s.basisType     = 'Continuous';  % Continuous / Discontinuous
s.coarseOrder   = 1;
s.dirac         = true;

%% COMPUTE
training = TrainingPreconditioning(s);
results  = training.compute();