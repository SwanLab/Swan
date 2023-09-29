classdef Filter_PDE_Density < handle
    
    properties (Access = private)
        mesh
        LHStype
        scale
        epsilon
        Anodal2Gauss
        problemLHS
        LHS
        bc
        quadrature
    end

    methods (Access = public)

        function obj = Filter_PDE_Density(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.computeBoundaryConditions();
            obj.epsilon = cParams.mesh.computeMeanCellSize();
            obj.Anodal2Gauss = obj.computeA();
            lhs = obj.createProblemLHS(cParams);
            obj.LHS = decomposition(lhs);
        end

        function x0 = getP0Function(obj,f,quadType)
            xRP1 =  obj.getP1fromP1(f,quadType);
            xR   = xRP1.fValues;
            x0   = zeros(obj.mesh.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*xR;
            end
        end

        function RHS = computeRHS(obj,f,quadType)
            s.fun = f;
            s.quadType = quadType;
            s.trial    = P1Function.create(obj.mesh, 1);
            s.type     = 'functionWithShapeFunction';
            s.mesh     = obj.mesh;
            in        = RHSintegrator.create(s);
            RHS          = in.RHS;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                lhs = obj.computeLHS();
                obj.LHS = decomposition(lhs);
            end
        end

        function xReg = getP1fromP1(obj,f,quadType)
            xReg = obj.getP1Function(f,quadType);
        end

        function xReg = getP1Function(obj,f,quadType)
            RHS       = obj.computeRHS(f,quadType);
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

        function A_nodal_2_gauss = computeA(obj)
            p1f = P1Function.create(obj.mesh,1);
            s.nnode   = obj.mesh.nnodeElem;
            s.nelem   = obj.mesh.nelem;
            s.npnod   = obj.mesh.nnodes;
            s.ngaus   = obj.quadrature.ngaus;
            s.connec  = obj.mesh.connec;
            s.shape   = p1f.computeShapeFunctions(obj.quadrature);
            Acomp = Anodal2gausComputer(s);
            Acomp.compute();
            A_nodal_2_gauss = Acomp.A_nodal_2_gauss;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
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

    end

end