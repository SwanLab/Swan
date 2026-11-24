function s = latticeSidesSettings()
s = latticeSettings();
s.hole.mU = 0.99;
for name = {'b','rho','t'}
    if isfield(s.dataset,name{1})
        s.dataset = rmfield(s.dataset,name{1});
    end
end
s.dataset.m          = s.hole.mU*(1-cos(linspace(0,pi,33)))/2;
s.dataset.a          = linspace(0,0.8,17);
s.dataset.symmetries = {'mirrorA','swapSides'};
s.dataset.folder     = 'latticeSidesDataset';
end