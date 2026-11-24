function s = latticeXiSettings()
s = latticeSettings();
s.hole.mU            = 0.99;
s.dataset            = rmfield(s.dataset,'b');
s.dataset.a          = linspace(-0.8,0.8,33);
s.dataset.t          = linspace(0,1,21);
s.dataset.symmetries = {};
s.dataset.folder     = 'latticeXiDataset';
end