function s = latticeSettings()
s.material       = struct('E',1,'nu',1/3,'voidRatio',1e-10);
s.mesh           = struct('divUnit',160,'strategy','Adapted');
s.transformation = struct('type','BoostSqueeze');
s.motif.type          = 'Square';
s.motif.pnorm         = 'Inf';
s.motif.center        = [0.5,0.5];
s.motif.tolerance     = 1e-4;
s.motif.maxIterations = 40;
s.motif.scaleBracket  = [1e-10,1];
s.dataset.rho        = linspace(0,1,20);
s.dataset.a          = tanh(linspace(0,atanh(0.98),30));
s.dataset.b          = exp(linspace(0,log(3),13));
s.dataset.symmetries = {'mirrorA','inverseB'};
s.dataset.nWorkers   = 4;
s.dataset.folder     = 'latticeDataset';
s.dataset.outputFile = 'dataset.mat';
end