classdef NewFemComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = NewFemComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            s = obj.createFEMparameters();
            obj.computation = NewFEM.create(s);
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
        
        function createMaterialProperties(obj)
            q = Quadrature.set(obj.computation.mesh.type);
            q.computeQuadrature('LINEAR');
            I = ones(obj.computation.mesh.nelem,q.ngaus);
            p.kappa = .9107*I;
            p.mu    = .3446*I;
            obj.computation.setMatProps(p);
        end

    end
end