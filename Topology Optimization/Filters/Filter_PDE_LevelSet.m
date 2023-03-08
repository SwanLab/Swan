classdef Filter_PDE_LevelSet < Filter

    properties (Access = private)
        levelSet
        Acomp
        Anodal2Gauss
        epsilon
        x_reg
        LHS
        bc
    end

    methods (Access = public)

        function obj = Filter_PDE_LevelSet(cParams)
            obj.init(cParams);
            obj.computeBoundaryConditions();
            obj.levelSet = cParams.designVariable;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
            obj.Anodal2Gauss = obj.computeA();
            lhs = obj.createProblemLHS();
            obj.LHS = decomposition(lhs);
        end

        function RHS = integrate_L2_function_with_shape_function(obj,x)
            ls = obj.levelSet.value;
            F = ones(size(ls));
            RHS = obj.computeRHS(F);
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
            x0 = zeros(obj.mesh.nelem,obj.field.quadrature.ngaus);
            for igaus = 1:obj.field.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*obj.x_reg;
            end
        end

        function x_reg = getP1fromP0(obj,x0)
            s.geometry = obj.field.geometry;
            s.x        = x0;
            RHS        = obj.Acomp.integrateP1FunctionWithShapeFunction(s);
            x_reg      = obj.solveFilter(RHS);
        end

    end

    methods (Access = private)

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
            s.ngaus   = obj.field.quadrature.ngaus;
            s.connec  = obj.mesh.connec;
            s.shape   = obj.field.interpolation.shape;
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
            s.field    = obj.field;
            s.type     = obj.LHStype;
            problemLHS = LHSintegrator.create(s);
            lhs        = problemLHS.compute(obj.epsilon);
            obj.bc     = obj.field.translateBoundaryConditions(obj.bc);
            lhs        = obj.bc.fullToReducedMatrix(lhs);
        end
        
        function fInt = computeRHSinBoundary(obj,fNodes)
            ls = obj.levelSet.value;
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                uMesh = obj.levelSet.getUnfittedMesh();
                s.mesh = uMesh;
                s.type = 'Unfitted';
                int = RHSintegrator.create(s);
                fInt = int.integrateInBoundary(fNodes);
            end
        end

    end

end