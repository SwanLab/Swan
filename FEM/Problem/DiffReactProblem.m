classdef DiffReactProblem < handle
    
    properties (GetAccess = public, SetAccess = protected)
        variables
    end
    
    properties (Access = private)
        dim
        mesh
        solver
        epsilon
        LHStype
        problemLHS
        problemData
        boundaryConditions
    end

    methods (Access = public)
        
        function obj = DiffReactProblem(cParams)
            obj.init(cParams);
            obj.computeDimensions();
            obj.createBoundaryConditions();
            obj.createSolver();
            obj.createProblemLHS();
        end
        
        function computeVariables(obj,rhs)
            bc  = obj.boundaryConditions;
            RHS = bc.fullToReducedVector(rhs);
            LHS = obj.computeLHS();
            x = obj.solver.solve(LHS,RHS);
            obj.variables.x = bc.reducedToFullVector(x);
        end
        
        function LHS = computeLHS(obj)
            lhs = obj.problemLHS.compute(obj.epsilon);
            LHS = obj.boundaryConditions.fullToReducedMatrix(lhs);
        end

        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
        end

        function M = getM(obj)
            M = obj.problemLHS.M;
        end

        function K = getK(obj)
            K = obj.problemLHS.K;
        end

        function dvol = computeDvolume(obj)
            int = obj.mesh.interpolation;
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            dvol = g.dvolu;
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.problemData.pdim = '1D';
            obj.problemData.scale = cParams.scale;
            obj.setLHStype(cParams);
        end

        function computeDimensions(obj)
            s.ngaus = [];
            s.mesh  = obj.mesh;
            s.pdim  = obj.problemData.pdim;
            dims    = DimensionVariables(s);
            dims.compute();
            obj.dim = dims;
        end

        function createBoundaryConditions(obj)
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.scale        = obj.problemData.scale;
            s.bc.dirichlet = [];
            s.bc.pointload = [];
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end
        
        function createSolver(obj)
            obj.solver = Solver.create();
        end

        function createProblemLHS(obj)
            s.type         = obj.LHStype;
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = [];
            obj.problemLHS = LHSintegrator.create(s);
        end

        function setLHStype(obj, cParams)
            isRobinAdded = isfield(cParams, 'isRobinTermAdded') ...
                               && cParams.isRobinTermAdded == 1;
            switch isRobinAdded
                case true
                    type = 'DiffReactRobin';
                case false
                    type = 'DiffReactNeumann';
            end
            obj.LHStype = type;
        end
    
    end

end