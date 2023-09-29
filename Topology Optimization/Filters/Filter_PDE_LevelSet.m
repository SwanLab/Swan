classdef Filter_PDE_LevelSet < handle

    properties (Access = private)
        mesh
        LHStype
        scale
        levelSet
        Anodal2Gauss
        epsilon
        problemLHS
        LHS
        bc
        quadrature
    end

    methods (Access = public)

        function obj = Filter_PDE_LevelSet(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.computeBoundaryConditions();
            obj.levelSet = cParams.designVariable;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
            obj.Anodal2Gauss = obj.computeA();
            lhs = obj.createProblemLHS(cParams);
            obj.LHS = decomposition(lhs);
        end

        function RHS = computeRHS(obj,f,quadType)
            % given design variable...
            %ls = obj.levelSet.value;
            ls = f.fValues;
            F = ones(size(ls));
            RHS = obj.computeRHSProjection(F);
        end

        function RHS = computeRHSoriginal(obj,f,quadType)
            s.quadType = quadType;
            s.fun      = f;
            s.trial    = P1Function.create(obj.mesh, 1);
            s.type     = 'functionWithShapeFunction';
            s.mesh     = obj.mesh;
            in        = RHSintegrator.create(s);
            RHS          = in.RHS;
        end

        function RHS = integrateFunctionAlongFacets(obj,F)
            RHS = obj.computeRHSinBoundary(F);
        end

        function x_reg = regularize(obj,F)
            RHS = obj.integrateFunctionAlongFacets(F);
            x_reg = obj.solveFilter(RHS);
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                lhs = obj.computeLHS();
                obj.LHS = decomposition(lhs);
            end
        end

        function xReg = getP1fromP1(obj,f,quadType)
            RHS       = obj.computeRHS(f,quadType);
            xR        = obj.solveFilter(RHS);
            p.fValues = xR;
            p.mesh    = obj.mesh;
            xReg      = P1Function(p);
        end

        function x0 = getP0Function(obj,f,quadType)
            xRP1 =  obj.getP1fromP1(f,quadType);
            xR   = xRP1.fValues;
            x0 = zeros(obj.mesh.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*xR;
            end
        end

        function xReg = getP1Function(obj,f,quadType)
            RHS       = obj.computeRHSoriginal(f,quadType);
            xR        = obj.solveFilter(RHS);
            p.fValues = xR;
            p.mesh    = obj.mesh;
            xReg      = P1Function(p);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.scale   = cParams.scale;
            if isfield(cParams,'LHStype')
                obj.LHStype = cParams.LHStype;
            else
                obj.LHStype = 'DiffReactNeumann';
            end
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end

        function fInt = computeRHSProjection(obj,fNodes)
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
            s.type = 'ShapeFunction';
            int = RHSintegrator.create(s);
        end

        function A = computeA(obj)
            p1f = P1Function.create(obj.mesh,1);
            s.nnode   = obj.mesh.nnodeElem;
            s.nelem   = obj.mesh.nelem;
            s.npnod   = obj.mesh.nnodes;
            s.ngaus   = obj.quadrature.ngaus;
            s.connec  = obj.mesh.connec;
            s.shape   = p1f.computeShapeFunctions(obj.quadrature);
            Acomp = Anodal2gausComputer(s);
            Acomp.compute();
            A = Acomp.A_nodal_2_gauss;
        end

        function x_reg = solveFilter(obj,RHS)
            RHS = obj.bc.fullToReducedVector(RHS);
            s.type = 'DIRECT';
            Solv = Solver.create(s);
            x = Solv.solve(obj.LHS,RHS);
            x_reg = obj.bc.reducedToFullVector(x);
        end

        function computeBoundaryConditions(obj)
            s.scale        = obj.scale;
            s.mesh         = obj.mesh;
            s.bc{1}.dirichlet = [];
            s.bc{1}.pointload = [];
            s.bc{1}.ndimf     = [];
            s.bc{1}.ndofs     = [];
            s.ndofs        = obj.mesh.nnodes;
            obj.bc         = BoundaryConditions(s);
            obj.bc.compute();
        end

        function lhs = createProblemLHS(obj,s)
            s.mesh         = obj.mesh;
            s.type         = obj.LHStype;
            obj.problemLHS = LHSintegrator.create(s);
            lhs            = obj.computeLHS();
        end

        function lhs = computeLHS(obj)
            lhs = obj.problemLHS.compute(obj.epsilon);
            lhs = obj.bc.fullToReducedMatrix(lhs);
        end

        function fInt = computeRHSinBoundary(obj,fNodes)
            ls = obj.levelSet.value;
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                uMesh = obj.levelSet.getUnfittedMesh();
                s.mesh = uMesh;
                s.type = 'ShapeFunction';
                int = RHSintegrator.create(s);
                fInt = int.integrateInBoundary(fNodes);
            end
        end

    end

end