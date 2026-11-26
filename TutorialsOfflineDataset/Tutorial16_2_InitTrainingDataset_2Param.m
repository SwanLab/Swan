% InitTrainingDataset
clc; clear;

%% COMMON INPUTS

% MESH
s.meshSource = 'Generated';      % Generated / Imported
s.ndim       = 2;
s.nelem      = [30 30];          % Just for generated mesh 
s.bounds     = [2 2];            % Just for generated mesh

% MATERIAL
s.inclusionType = 'Weak';        % Void /Weak
s.young         = 1;
s.poisson       = 1/3;

% TRAINING
s.method        = 'Multiscale';  % Multiscale / Eifem
s.nSubdomains   = [5 5];
s.domainIndices = [3 3];
s.basisType     = 'Continuous';  % Continuous / Discontinuous
s.coarseOrder   = 1;
s.dirac         = true;

%% GEOMETRY
g.type        = 'CrossedSquare';
g.length      = 2;
g.xCoorCenter = 1;
g.yCoorCenter = 1;

%% PARAMETERS
tFrameValues = 0.1:0.05:0.45;
tCrossValues = 0.1:0.05:0.6;

nSamples = numel(tFrameValues)*numel(tCrossValues);

cases           = cell(nSamples, 1);
parameterValues = zeros(nSamples, 2);

%% DEFINE SAMPLES
count = 0;

for i = 1:numel(tFrameValues)
    for j = 1:numel(tCrossValues)
        count = count + 1;
        g.tFrame = tFrameValues(i);
        g.tCross = tCrossValues(j);
        sampleInputs = s;
        sampleInputs.geometry = GeometricalFunction(g);

        % Store inputs and matching parameter values
        cases{count} = sampleInputs;
        parameterValues(count,:) = [g.tFrame, g.tCross];
    end
end

%% DATASET INPUTS
p.cases           = cases;
p.parameterNames  = {'tFrame', 'tCross'};
p.parameterValues = parameterValues;

%% EXPORT OPTIONS
p.directory = fullfile('TutorialsOfflineDataset/Results', 'Lattice', s.method, s.inclusionType);
p.saveMAT   = true;
p.saveCSV   = true;

%% COMPUTE AND EXPORT
dataset = TrainingDataset(p);
samples = dataset.compute();
dataset.export(samples);