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
            s = FemInputReader_GiD().read(obj.testName);
            obj.computation = NewFEM.create(s);
            obj.computation.computeVariables();
        end
    end

    methods (Access = private)

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