classdef NewFilter_PDE_LevelSet < handle

    properties (Access = private)
        interp
        levelSet
        Acomp
        Anodal2Gauss
        epsilon
        x_reg
        quadrature
        geometry
        dim
        LHS
        bc
        field
    end

    properties (Access = private)
        mesh
        quadratureOrder
        femSettings
        LHStype
    end

    methods (Access = public)

        function obj = NewFilter_PDE_LevelSet(cParams)
            obj.init(cParams);
            obj.computeDimension();
            obj.computeBoundaryConditions();
            obj.levelSet = cParams.designVariable;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function preProcess(obj)
            s.mesh            = obj.mesh;
            s.quadratureOrder = obj.quadratureOrder;
            P1proc            = P1preProcessor(s);
            P1proc.preProcess();
            obj.storeParams(P1proc);
            obj.Anodal2Gauss = obj.computeA();
            lhs = obj.createProblemLHS();
            obj.LHS = decomposition(lhs);
        end

        function RHS = integrate_L2_function_with_shape_function(obj,x)
            ls = obj.levelSet.value;
            F = ones(size(ls));
            RHS = obj.computeRHS(F);
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                lhs = obj.createProblemLHS();
                obj.LHS = decomposition(lhs);
            end
        end

        function x_reg = getP1fromP1(obj,x)
            RHS = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solveFilter(RHS);
        end

        function x0 = getP0fromP1(obj,x)
            obj.x_reg =  obj.getP1fromP1(x);
            x0 = zeros(obj.mesh.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*obj.x_reg;
            end
        end

        function x_reg = getP1fromP0(obj,x0)
            s.geometry = obj.geometry;
            s.x        = x0;
            RHS        = obj.Acomp.integrateP1FunctionWithShapeFunction(s);
            x_reg      = obj.solveFilter(RHS);
        end

        function q = getQuad(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
            obj.femSettings = cParams.femSettings;
            if isfield(cParams.femSettings,'LHStype')
                obj.LHStype = cParams.femSettings.LHStype;
            else
                obj.LHStype = 'DiffReactNeumann';
            end
            obj.createField();
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.field = Field(s);
        end

        function storeParams(obj,P1proc)
            obj.quadrature = P1proc.quadrature;
            obj.interp     = P1proc.interp;
            obj.geometry   = P1proc.geometry;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end

        function fInt = computeRHS(obj,fNodes)
            ls = obj.levelSet.value;
            int  = obj.obtainRHSintegrator();
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                fInt = int.integrateInDomain(fNodes);
            end
        end

        function int = obtainRHSintegrator(obj)
            uMesh = obj.levelSet.getUnfittedMesh();
            s.mesh = uMesh;
            s.type = 'Unfitted';
            int = RHSintegrator.create(s);
        end

        function A = computeA(obj)
            s.nnode   = obj.mesh.nnodeElem;
            s.nelem   = obj.mesh.nelem;
            s.npnod   = obj.mesh.nnodes;
            s.ngaus   = obj.quadrature.ngaus;
            s.connec  = obj.mesh.connec;
            s.shape   = obj.interp.shape;
            obj.Acomp = Anodal2gausComputer(s);
            obj.Acomp.compute();
            A = obj.Acomp.A_nodal_2_gauss;
        end

        function x_reg = solveFilter(obj,RHS)
            RHS = obj.bc.fullToReducedVector(RHS);
            s.type = 'DIRECT';
            Solv = Solver.create(s);
            x = Solv.solve(obj.LHS,RHS);
            x_reg = obj.bc.reducedToFullVector(x);
        end

        function computeDimension(obj)
            s.name  = 'x';
            s.mesh  = obj.mesh;
            obj.dim = DimensionScalar(s);
        end

        function computeBoundaryConditions(obj)
            s.dim          = obj.dim;
            s.scale        = obj.femSettings.scale;
            s.mesh         = obj.mesh;
            s.bc{1}.dirichlet = [];
            s.bc{1}.pointload = [];
            s.bc{1}.ndimf     = [];
            s.bc{1}.ndofs     = [];
            s.ndofs        = obj.dim.ndofs;
            obj.bc         = BoundaryConditions(s);
            obj.bc.compute();
        end

        function lhs = createProblemLHS(obj)
            s      = obj.femSettings;
            s.type = obj.LHStype;
            s.dim  = obj.dim;
            s.mesh = obj.mesh;
            s.globalConnec = [];
            problemLHS = LHSintegrator.create(s);
            lhs = problemLHS.compute(obj.epsilon);
            lhs = obj.bc.fullToReducedMatrix(lhs);
        end

    end

end