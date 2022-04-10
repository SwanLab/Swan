%% Hyperelastic
% filename = 'test_hyperelastic';
% s = createFEMparameters(filename);
% 
% fem = FEM.create(s);
%% Thermal
filename = 'Cantilever119Kelem';
s = createFEMparameters(filename);
fem = FEM.create(s);
tic
fem.solve();
toc
fem.print(filename)


%% Thermal
% filename = 'test_thermal';
% s = createFEMparameters(filename);
% fem = FEM.create(s);
% fem.solve();
% fem.print(filename)

%% Functions
function s = createFEMparameters(file)
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