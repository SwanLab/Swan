classdef StructuralTrussProblem < handle
    
    properties (Access = public)
        stress
        displacement
    end

    properties (Access = private) % In
        coord
        connec
        material
    end

    properties (Access = private) % Calc
        solver
        LHS, RHS
    end

    methods (Access = public)

        function obj = StructuralTrussProblem(cParams)
            obj.init(cParams);
            obj.createSolver();
        end

        function solve(obj)
            obj.computeLHS();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.coord    = cParams.coord;
            obj.connec   = cParams.connec;
            obj.material = cParams.material;
        end

        function createSolver(obj)
            s.type = 'DIRECT';
            obj.solver = Solver.create(s);
        end

        function computeLHS(obj)
            % To optimize by calculating lengths just once
            % and just multiply the LHS by E*A/L
            s.coord    = obj.coord;
            s.connec   = obj.connec;
            s.material = obj.material;
            lhs = LHSintegrator_StiffnessBeam(s);
            lhs.compute();
        end

        function computeRHS(obj)
        end

    end

end
