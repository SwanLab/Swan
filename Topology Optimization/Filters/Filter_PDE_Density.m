classdef Filter_PDE_Density < Filter
    
    properties (Access = private)
        epsilon
        Anodal2Gauss
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
            lhs = obj.createProblemLHS();
            obj.LHS = decomposition(lhs);
        end

        function x0 = getP0fromP1(obj,x)
            xR =  obj.getP1fromP1(x);
            x0 = zeros(obj.mesh.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*xR;
            end
        end

        function RHS = integrate_L2_function_with_shape_function(obj,x)
            p.fValues = x;
            p.mesh = obj.mesh;
            f = P1Function(p);
            s.fun = f;
            s.quadType = 'QUADRATICMASS';
            RHS = obj.computeRHS(s);
        end

        function RHS = computeRHS(obj,cParams)
            s.quadType = cParams.quadType;
            s.fun      = cParams.fun;
            s.trial    = P1Function.create(obj.mesh, 1);
            s.type     = 'functionWithShapeFunction';
            s.mesh     = obj.mesh;
            in        = RHSintegrator.create(s);
            RHS          = in.RHS;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                lhs = obj.createProblemLHS();
                obj.LHS = decomposition(lhs);
            end
        end

        function xReg = getP1fromP1(obj,x)
            RHS  = obj.integrate_L2_function_with_shape_function(x);
            xReg = obj.solveFilter(RHS);
        end

        function xReg = getP1fromP0(obj,x0)
            nelem     = size(x0,1);
            ngaus     = size(x0,2);
            s.fValues = reshape(x0',[1,ngaus,nelem]);
            s.mesh    = obj.mesh;
            s.quadrature = obj.quadrature;
            f         = FGaussDiscontinuousFunction(s);
            a.quadType = 'LINEAR'; % may be an input !!
            a.fun = f;
            RHS = obj.computeRHS(a);
            xReg      = obj.solveFilter(RHS);
        end

    end

    methods (Access = private)

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
            s.scale        = obj.femSettings.scale;
            s.mesh         = obj.mesh;
            s.bc{1}.dirichlet = [];
            s.bc{1}.pointload = [];
            s.bc{1}.ndimf     = [];
            s.bc{1}.ndofs     = [];
            s.ndofs        = obj.mesh.nnodes;
            obj.bc         = BoundaryConditions(s);
            obj.bc.compute();
        end

        function lhs = createProblemLHS(obj)
            s          = obj.femSettings;
            s.mesh     = obj.mesh;
            s.type     = obj.LHStype;
            problemLHS = LHSintegrator.create(s);
            lhs        = problemLHS.compute(obj.epsilon);
            lhs        = obj.bc.fullToReducedMatrix(lhs);
        end

    end

end