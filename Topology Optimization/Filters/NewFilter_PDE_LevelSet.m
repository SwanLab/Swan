classdef NewFilter_PDE_LevelSet < handle

    properties (Access = public)
        unfittedMesh
        diffReacProb
        ngaus
        nelem
    end

    properties(Access = private)
        integrator
        domainType
        interp
        levelSet
        Anodal2Gauss
        LHS
        epsilon
        x
        x_reg
        mesh
        quadratureOrder
        npnod
        quadrature
        geometry
        nnode
        shape
    end

    methods (Access = public)

        function obj = NewFilter_PDE_LevelSet(cParams)
            obj.init(cParams)
            obj.levelSet = cParams.designVariable;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
            obj.domainType = cParams.domainType;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.computeGeometry();
            obj.nelem = obj.mesh.nelem;
            obj.npnod = obj.mesh.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.Anodal2Gauss = obj.computeA();
        end

        function preProcess(obj)
            s.mesh            = obj.mesh;
            s.quadratureOrder = obj.quadratureOrder;
            P1proc            = P1preProcessor(s);
            P1proc.preProcess();
            obj.storeParams(P1proc);
            obj.Anodal2Gauss = obj.computeA();
            obj.diffReacProb.setEpsilon(obj.epsilon);
            obj.computeLHS();
        end

        function RHS = integrate_L2_function_with_shape_function(obj,x)
            ls = obj.levelSet.value;
            F = ones(size(ls));
            RHS = obj.computeRHS(F);
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                obj.diffReacProb.setEpsilon(epsilon);
                obj.computeLHS();
            end
        end

        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end

        function x_reg = getP1fromP1(obj,x)
            RHS = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solve_filter(RHS);
        end

        function RHS = integrate_function_along_facets(obj,F)
            RHS = obj.computeRHSinBoundary(F);
        end

        function obj = createDiffReacProblem(obj,cParams)
            s = cParams.femSettings;
            if isprop(cParams,'mesh')
                s.mesh = cParams.mesh;
            end
            switch s.scale
                case 'MACRO'
                    obj.diffReacProb = NewDiffReactProblem(s);
                case 'MICRO'
                    obj.diffReacProb = NewDiffReactProblemMicro(s);
            end
        end

        function x0 = getP0fromP1(obj,x)
            obj.x_reg =  obj.getP1fromP1(x);
            x0 = zeros(obj.mesh.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*obj.x_reg;
            end
        end

        function x_reg = getP1fromP0(obj,x0)
            RHS = obj.integrate_P1_function_with_shape_function(x0);
            x_reg = obj.solve_filter(RHS);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.createDiffReacProblem(cParams);
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function storeParams(obj,P1proc)
            obj.quadrature = P1proc.quadrature;
            obj.interp     = P1proc.interp;
            obj.geometry   = P1proc.geometry;
            obj.nnode      = obj.mesh.nnode;
            obj.shape      = obj.interp.shape;
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

        function fInt = computeRHSinBoundary(obj,fNodes)
            ls = obj.levelSet.value;
            int  = obj.obtainRHSintegrator();
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                fInt = int.integrateInBoundary(fNodes);
            end
        end

        function int = obtainRHSintegrator(obj)
            uMesh = obj.levelSet.getUnfittedMesh();
            s.mesh = uMesh;
            s.type = 'Unfitted';
            int = Integrator.create(s);
        end

        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.type);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end

        function createInterpolation(obj)
            obj.interp = Interpolation.create(obj.mesh,'LINEAR');
        end

        function computeGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interp);
        end

        function A_nodal_2_gauss = computeA(obj)
            s.nnode  = obj.nnode;
            s.nelem  = obj.nelem;
            s.npnod  = obj.npnod;
            s.ngaus  = obj.ngaus;
            s.connec = obj.mesh.connec;
            s.shape  = obj.shape;
            Acomp = Anodal2gausComputer(s);
            Acomp.compute();
            A_nodal_2_gauss = Acomp.A_nodal_2_gauss;
        end

        function computeLHS(obj)
            lhs = obj.diffReacProb.computeLHS();
            obj.LHS = decomposition(lhs);
        end

        function x_reg = solve_filter(obj,RHS)
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end

        function intX = integrate_P1_function_with_shape_function(obj,x)
            ndof = size(obj.Anodal2Gauss{1},2);
            intX = zeros(ndof,1);
            for igaus = 1:obj.quadrature.ngaus
                dVG = obj.geometry.dvolu(:,igaus);
                xG = x(:,igaus);
                A = obj.Anodal2Gauss{igaus};
                intX = intX + A'*(xG.*dVG);
            end
        end

    end

end