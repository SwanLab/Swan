classdef MicroComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = MicroComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj) %new
            s = obj.createFEMparameters();
            femSolver = NewElasticProblemMicro(s);
            femSolver.solve();
            femSolver.computeChomog();
            obj.computation = femSolver;
        end

%         function compute(obj) %old
%             femSolver = Elastic_Problem_Micro.create(obj.testName);
%             props.kappa = .9107;
%             props.mu    = .3446;
%             femSolver.setMatProps(props);
%             femSolver.computeChomog();
%             obj.computation = femSolver;
%         end
    end

    methods (Access = private)

        function s = createFEMparameters(obj)
            gidParams = obj.createGiDparameters();
            s.dim       = gidParams.pdim;
            s.type      = gidParams.ptype;
            s.scale     = gidParams.scale;
            s.mesh      = gidParams.mesh;
            s.dirichlet = gidParams.dirichlet;
            s.pointload = gidParams.pointload;
            if isequal(s.scale,'MICRO')
                s.masterSlave = gidParams.masterSlave;
            end
        end

        function gidParams = createGiDparameters(obj)
            file = obj.testName;
            gidReader = FemInputReader_GiD();
            gidParams = gidReader.read(file);
        end
        
    end

end

