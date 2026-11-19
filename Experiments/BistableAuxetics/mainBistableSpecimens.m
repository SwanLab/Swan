if ~exist('refinement','var')
    refinement = 1;
end
if ~exist('strainStep','var')
    strainStep = 0.005;
end
if ~exist('contact','var')
    contact = 'barrier';
end
if ~exist('regularization','var')
    regularization = 1e-3;
end
if ~exist('mediumStiffness','var')
    mediumStiffness = 1e-6;
end

here = fileparts(mfilename('fullpath'));

for iRun = 1:size(runs,1)
    s.geometry        = createGeometry(runs{iRun,1});
    s.name            = createRunName(runs{iRun,1},runs{iRun,2},refinement,contact,regularization,mediumStiffness);
    s.refinement      = refinement;
    s.conditionsType  = runs{iRun,2};
    s.strains         = 0:strainStep:runs{iRun,3};
    s.matProp         = createNeohookeanProperties(1,0.45);
    s.tolerance       = 1e-10;
    s.maxIter         = 100;
    s.contact         = contact;
    s.regularization  = regularization;
    s.mediumStiffness = mediumStiffness;
    s.outDir          = createOutputFolder(here,s.name);

    t0 = tic;
    computer = BistableSpecimenComputer(s);
    computer.plotSetup();
    computer.compute();
    r = computer.results;
    r.wallTime = toc(t0);
    save(fullfile(s.outDir,[s.name,'.mat']),'r');
    printSummary(r,s);
end

function g = createGeometry(specimen)
    nCells = str2double(specimen(end));
    switch specimen(1:end-1)
        case 'hexagon'
            g = PinwheelSpecimenGeometry(PinwheelSpecimenGeometry.hexagonParams(nCells));
        case 'hexagonPrint'
            g = PinwheelSpecimenGeometry(PinwheelSpecimenGeometry.hexagonPrintParams(nCells));
        case 'squareThick'
            g = PinwheelSpecimenGeometry(PinwheelSpecimenGeometry.squareParams(nCells,0.05));
        case 'squareThin'
            g = PinwheelSpecimenGeometry(PinwheelSpecimenGeometry.squareParams(nCells,0.025));
        case 'curvedBeam'
            g = ImportedSpecimenGeometry(ImportedSpecimenGeometry.curvedBeamParams(nCells));
        case 'frame'
            g = FrameSpecimenGeometry(FrameSpecimenGeometry.paperFrameParams());
    end
end

function name = createRunName(specimen,type,refinement,contact,regularization,mediumStiffness)
    name = sprintf('%s_%s_r%g',specimen,type,refinement);
    if ~strcmp(contact,'none')
        name = [name,'_',contact];
    end
    if strcmp(contact,'thirdMedium')
        name = sprintf('%s_kr%g_g%g',name,regularization,mediumStiffness);
    end
end

function p = createNeohookeanProperties(mu,nu)
    lambda3D = 2*mu*nu/(1-2*nu);
    p.mu     = mu;
    p.lambda = 2*lambda3D*mu/(lambda3D+2*mu);
end

function outDir = createOutputFolder(here,name)
    outDir = fullfile(here,'results',name);
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
end

function printSummary(r,s)
    fprintf('\n=== %s: nodes=%d elems=%d time=%.0fs lastStrain=%.3f of %.3f\n',s.name, ...
        r.nNodes,r.nElem,r.wallTime,r.strain(end),s.strains(end));
    disp(table(r.strain,r.energy,r.force,r.residual,r.minJ,r.maxStretch,r.overlaps,r.iter, ...
        'VariableNames',{'strain','energy','force','residual','minJ','maxStretch','overlaps','iter'}));
end
