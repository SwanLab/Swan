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
            obj.createProblemLHS(cParams);
        end

        function computeVariables(obj,rhs)
            bc  = obj.boundaryConditions;
            RHS = bc.fullToReducedVector(rhs);
            LHS = obj.computeLHS(obj.epsilon);
            x = obj.solver.solve(LHS,RHS);
            obj.variables.x = bc.reducedToFullVector(x);
        end
        
        function LHS = computeLHS(obj, epsilon)
            obj.epsilon = epsilon;
            lhs = obj.problemLHS.compute(epsilon);
            LHS = obj.boundaryConditions.fullToReducedMatrix(lhs);
        end
       
        function print(obj,filename)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            s.quad = quad;
            s.mesh = obj.mesh;
            s.iter = 0;
            s.fields    = obj.variables.x;
            s.ptype     = 'DIFF-REACT';
            s.ndim      = 3;
            s.pdim      = obj.problemData.pdim;
            s.type      = 'ScalarNodal';
            fPrinter = FemPrinter(s);
            fPrinter.print(filename);
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh              = cParams.mesh;
            obj.LHStype           = cParams.LHStype;
            obj.problemData.pdim  = '1D';
            obj.problemData.scale = cParams.scale;
        end

        function computeDimensions(obj)
            s.name = 'x';
            s.mesh = obj.mesh;
            dims   = DimensionScalar(s);
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

        function createProblemLHS(obj,cParams)
            s.type         = obj.LHStype;
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            if isfield(cParams,'isAnisotropyAdded')
                s.isAnisotropyAdded = cParams.isAnisotropyAdded;
            end
            if isfield(cParams,'CAnisotropic')
                s.CAnisotropic = cParams.CAnisotropic;
            end
            s.globalConnec = [];
            obj.problemLHS = LHSintegrator.create(s);
        end
    
    end

end