function DehomogenisingCases


%iteration = 20;
%fCase = 'ExperimentingPlot';
%fCase = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';
%folder = ['/home/alex/git-repos/Swan/Output/',fCase];

% iteration = 235;
% fCase = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';
% folder = '/media/alex/My Passport/LatticeResults/CantileverSymmetricMeshSuperEllipsePDEGradientEpsilonH';


%iteration = 312;
%fCase = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';
%folder = '/media/alex/My Passport/LatticeResults/CantileverSymmetricMeshSuperEllipsePDEGradientDensityEpsilonH';

%iteration = 250;
%fCase = 'ExperimentingPlot';
%folder = '/media/alex/My Passport/LatticeResults/StressNormRectangleRotationSmall';

%iteration = 296;
%fCase = 'ExperimentingPlot';
%folder = '/media/alex/My Passport/LatticeResults/StressNormSuperEllipseRotationSmall';

% iteration = 262;
% fCase = 'CantileverSymmetricFixingMaxStressZone';
% folder = '/media/alex/My Passport/LatticeResults/CantileverSymmetricCoarse/CantileverSymmetricFixingMaxStressZone';

folder = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
fCase = 'LshapeCoarseSuperEllipseDesignVariable';
iteration = 665;

% folder = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/ExampleCompliance/';  
% fCase = 'ExperimentingPlotSuperEllipse';
% iteration = 64;

s.fileName = [fCase,num2str(iteration)];
s.folderPath = fullfile(folder);

d = DehomogenisationPrinter(s);
d.print();
end

