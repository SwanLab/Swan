classdef DiffReactProblem < handle
    
    properties (GetAccess = public, SetAccess = protected)
        variables
    end

    % Filter_PDE_Iso_Total
    % e^2 IsoLaplac(u) + u = u* (neuman)    
    
    % Filter_PDE_Iso_Relative
    % e^2 IsoLaplac(u) + u = u* (Robin)
    
    % Filter_PDE_Ani_Total   
    % e^2 AniLaplac(u) + u = u* (neuman)

    % Filter_PDE_Ani_Total    
    % e^2 AniLaplac(u) + u = u* (Robin)

    % IsoLaplac(u) = f
    % IsoVectLaplac(u) = f
    
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
            obj.computeStiffnessMatrix(cParams);
            obj.computeMassMatrix();
            obj.createProblemLHS();
        end
        
        function solve(obj)
            obj.epsilon = .1857;
            bc  = obj.boundaryConditions;
            M = obj.problemLHS.K;
            rhs = M*ones(size(M,1), 1);
            RHS = bc.fullToReducedVector(rhs);
            LHS = obj.computeLHS();
            x = obj.solver.solve(LHS,RHS);
            obj.variables.x = bc.reducedToFullVector(x);
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
       
        function print(obj,filename)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            s.quad = quad;
            s.mesh = obj.mesh;
            s.iter = 0;
            s.fields    = obj.variables.x;
            s.ptype     = 'DIFF-REACT';
            s.ndim      = 2;
            s.pdim      = obj.problemData.pdim;
            s.type      = 'ScalarNodal';
            fPrinter = FemPrinter(s);
            fPrinter.print(filename);
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

        function computeStiffnessMatrix(obj,cParams)
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
        end

        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
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