% InitTrainingDataset
clc; clear;


%% COMMON INPUTS

% MESH
s.meshSource = 'Generated';       % Generated / Imported
s.ndim       = 2;
s.nelem      = [20 20];
s.bounds     = [2 2];

% MATERIAL
s.inclusionType = 'Weak';         % Void /Weak
s.young         = 1;
s.poisson       = 1/3;

% TRAINING
s.method        = 'Multiscale';   % Multiscale / Eifem
s.nSubdomains   = [5 5];
s.domainIndices = [3 3];
s.basisType     = 'Continuous';   % Continuous / Discontinuous
s.coarseOrder   = 1;
s.dirac         = true;

%% GEOMETRY
g.type        = 'CircleInclusion';
g.xCoorCenter = 1;
g.yCoorCenter = 1;

%% DEFINE SAMPLES
radii = 0.1:0.1:0.8;
cases = cell(numel(radii), 1);

for i = 1:numel(radii)
    g.radius = radii(i);

    sampleInputs = s;
    sampleInputs.geometry = GeometricalFunction(g);

    cases{i} = sampleInputs;
end

%% DATASET INPUTS
p.cases           = cases;
p.parameterNames  = {'radius'};
p.parameterValues = radii(:);

%% EXPORT OPTIONS
p.directory = fullfile('TutorialsOfflineDataset/Results', 'Circle', s.method, s.inclusionType);
p.saveMAT   = true;
p.saveCSV   = true;

%% COMPUTE AND EXPORT
dataset = TrainingDataset(p);
samples = dataset.compute();
dataset.export(samples);