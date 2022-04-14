%% Hyperelastic
% Microstructure
filename = 'test2d_micro';
%filename = 'IrrHexagon50x25x50';
s = createParameters(filename);
d = createDensity(s.mesh);


fem = FEM.create(s);
fem.computeChomog();
fem.print(filename);


function d = createDensity(mesh)
s.type = 'Density';
s.mesh = mesh;
s.value = [];
s.initialCase = 'circleInclusion';
s.creatorSettings.type = 'FromLevelSet';
s.creatorSettings.fracRadius = 0.5;
d = Density(s);

s.mesh = mesh;
s.field = d.value;
p = NodalFieldPlotter(s);
p.plot()
end


% %% Thermal
% filename = 'test_thermal';
% s = createFEMparameters(filename);
% fem = FEM.create(s);
% fem.solve();
% fem.print(filename)

%% Functions
function s = createParameters(file)
gidParams = createGiDparameters(file);
s.dim       = gidParams.pdim;
s.type      = gidParams.ptype;
s.scale     = gidParams.scale;
s.mesh      = gidParams.mesh;
s.dirichlet = gidParams.dirichlet;
s.pointload = gidParams.pointload;
end

function gidParams = createGiDparameters(file)
gidReader = FemInputReader_GiD();
gidParams = gidReader.read(file);
end

function postProcess(fileName)
dI = obj.createPostProcessDataBase(fileName);
postprocess = Postprocess('ScalarNodal',dI);
q = obj.element.quadrature;
d.fields = obj.variables;
d.quad = q;
postprocess.print(obj.iter,d);
end

function d = createPostProcessDataBase(obj,fileName)
dI.mesh    = obj.mesh;
dI.outName = fileName;
ps = PostProcessDataBaseCreator(dI);
d = ps.getValue();
end
