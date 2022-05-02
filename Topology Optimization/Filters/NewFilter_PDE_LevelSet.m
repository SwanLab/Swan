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
        LHS
    end

    properties (Access = private)
        mesh
        quadratureOrder
        femSettings
        LHStype
    end

    methods (Access = public)

        function obj = NewFilter_PDE_LevelSet(cParams)
            obj.init(cParams)
            obj.init2(cParams);
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

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
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
            int = Integrator.create(s);
        end

        function A = computeA(obj)
            s.nnode   = obj.mesh.nnode;
            s.nelem   = obj.mesh.nelem;
            s.npnod   = obj.mesh.npnod;
            s.ngaus   = obj.quadrature.ngaus;
            s.connec  = obj.mesh.connec;
            s.shape   = obj.interp.shape;
            obj.Acomp = Anodal2gausComputer(s);
            obj.Acomp.compute();
            A = obj.Acomp.A_nodal_2_gauss;
        end

        function x_reg = solveFilter(obj,RHS)
            x_reg = obj.LHS\(RHS);
        end
        % inicio ida de olla:
        function lhs = createProblemLHS(obj)
            s.type         = obj.LHStype;
            ss.name = 'x';
            ss.mesh = obj.mesh;
            dims   = DimensionScalar(ss);
            s.dim          = dims;
            s.mesh         = obj.mesh;
            if isfield(obj.femSettings,'isAnisotropyAdded')
                s.isAnisotropyAdded = obj.femSettings.isAnisotropyAdded;
            end
            if isfield(obj.femSettings,'CAnisotropic')
                s.CAnisotropic = obj.femSettings.CAnisotropic;
            end
            s.globalConnec = [];
            problemLHS = LHSintegrator.create(s);
            lhs = problemLHS.compute(obj.epsilon);
            s.scale = obj.femSettings.scale;
            s.bc.dirichlet = [];
            s.bc.pointload = [];
            bc = BoundaryConditions(s);
            bc.compute();
            lhs = bc.fullToReducedMatrix(lhs);
        end

        function init2(obj,cParams)
            obj.femSettings = cParams.femSettings;
            if isfield(cParams.femSettings,'LHStype')
                obj.LHStype = cParams.femSettings.LHStype;
            else
                obj.LHStype = 'DiffReactNeumann';
            end
        end

    end

end