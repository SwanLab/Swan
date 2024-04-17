classdef HyperelasticProblem < handle
    
    properties (Access = public)
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        quadrature
        boundaryConditions, BCApplier

        stiffness
        forces
        solver, solverType, solverMode, solverCase
        scale
        
        strain, stress
    end

    properties (Access = protected)
        mesh 
        material  
        displacementFun
    end

    methods (Access = public)

        function obj = HyperelasticProblem()
            obj.init();
            obj.createDisplacementFun();
            obj.computeInternalForces();
        end

        function solve(obj)
        end

    end

    methods (Access = private)

        function init(obj)
            obj.mesh = UnitHexaMesh(5,5,5);
            obj.material.lambda = 1;
            obj.material.mu = 1;

        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function computeInternalForces(obj)
            s.mesh = obj.mesh;
            test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            rhs = RHSintegrator_Hyperelasticity(s);
            rhs.compute(obj.uFun, test);
        end

    end

end
