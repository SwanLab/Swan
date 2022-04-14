%% Hyperelastic
% filename = 'test_hyperelastic';
% s = createFEMparameters(filename);
% 
% fem = FEM.create(s);
%% Performance
tests = {'Cantilever1Kelem', 'Cantilever4Kelem', 'Cantilever18Kelem', ...
    'Cantilever50Kelem', 'Cantilever74Kelem', 'Cantilever119Kelem', ...
    'Cantilever216Kelem', 'Cantilever338K'};
index = 1;
for test = tests
    s = createFEMparameters(test{1});
    fem = FEM.create(s);
    tic
    for i  = 1:5
        tic
        fem.solve();
        resultsNew{index}(i) = toc;
    end
    disp(resultsNew)
    index = index + 1;
end
a = 0;

% filename = 'Cantilever119Kelem';
% s = createFEMparameters(filename);
% fem = FEM.create(s);
% tic
% fem.solve();
% toc
% % fem.print(filename)


% Microstructure
%filename = 'test2d_micro';
filename = 'IrrHexagon50x25x50';
s = createFEMparameters(filename);
fem = FEM.create(s);
fem.computeChomog();
fem.print(filename);


% %% Thermal
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
s.masterSlave = gidParams.masterSlave;
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