classdef ThermalComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = ThermalComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            s = obj.createFEMparameters();
            obj.computation = FEM.create(s);
            obj.computation.solve();
        end
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
        end

        function gidParams = createGiDparameters(obj)
            file = obj.testName;
            gidReader = FemInputReader_GiD();
            gidParams = gidReader.read(file);
        end

    end
end