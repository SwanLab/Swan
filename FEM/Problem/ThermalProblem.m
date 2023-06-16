classdef ThermalProblem < handle
    
    properties (Access = public)
        temperature
    end

    properties (Access = private)
        boundaryConditions
        solver        
        LHS
        RHS
    end

    properties (Access = private)
        mesh        
        scale
        inputBC  
        rhsFun
    end

    methods (Access = public)

        function obj = ThermalProblem(cParams)
            obj.init(cParams);
            obj.createTemperatureFun();
            obj.createBoundaryConditions();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeLHS();
            obj.computeRHS();
            obj.computeTemperatures();
        end

        function print(obj,filename)
            [fun, funNames] = obj.getFunsToPlot();
            a.mesh     = obj.mesh;
            a.filename = filename;
            a.fun      = fun;
            a.funNames = funNames;
            pst = FunctionPrinter_Paraview(a);
            pst.print();
        end
        function [fun, funNames] = getFunsToPlot(obj)
            fun = {obj.temperature};
            funNames = {'Temperature'};
        end        

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh        = cParams.mesh;
            obj.scale       = cParams.scale;
            obj.inputBC     = cParams.bc;
            obj.rhsFun      = cParams.rhsFun;
        end

        function createTemperatureFun(obj)
            nDimf  = 1;
            t = P1Function.create(obj.mesh, nDimf);
            obj.temperature = t;
        end

        function createBoundaryConditions(obj)
            dim = obj.getFunDims();
            bc = obj.inputBC;
            bc.ndimf = dim.ndimf;
            bc.ndofs = dim.ndofs;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            s.bc    = {bc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end  

        function dim = getFunDims(obj)
            d.ndimf  = obj.temperature.ndimf;
            d.nnodes = size(obj.temperature.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createSolver(obj)
            s.type = 'DIRECT';
            obj.solver = Solver.create(s);
        end        

        function computeLHS(obj)
            s.type     = 'StiffnessMatrix';
            s.mesh     = obj.mesh;
            s.test     = obj.temperature;
            s.trial    = obj.temperature;
            lhs = LHSintegrator.create(s);
            obj.LHS = lhs.compute();
        end

        function computeRHS(obj)
            rhsV    = obj.computeVolumetricRHS();
            obj.RHS = rhsV;
        end

        function rhs = computeVolumetricRHS(obj)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.test     = obj.temperature ;
            RHSint     = RHSintegrator.create(s);
            f          = obj.rhsFun;
            rhs = RHSint.compute(f);            
        end

        function t = computeTemperatures(obj)
            bc = obj.boundaryConditions;
            Kred = bc.fullToReducedMatrix(obj.LHS);
            Fred = bc.fullToReducedVector(obj.RHS);
            t = obj.solver.solve(Kred,Fred);
            t = bc.reducedToFullVector(t);
            obj.temperature.fValues = t;
        end

    end

end