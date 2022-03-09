classdef FemComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = FemComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            obj.computation = FEM.create(obj.testName);
            obj.createMaterialProperties();
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
            obj.computation.setMatProps(p)
        end

    end
end