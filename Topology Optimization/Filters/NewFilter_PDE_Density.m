classdef NewFilter_PDE_Density < handle
    
    properties (Access = private)
        epsilon
        diffReacProb
        mesh
        quadratureOrder
        Anodal2Gauss
        quadrature
        M
        interp
        geometry
        x_reg
        LHS
    end

    methods (Access = public)

        function obj = NewFilter_PDE_Density(cParams)
            obj.init(cParams);
            obj.createMassMatrix();
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function preProcess(obj)
            s.mesh            = obj.mesh;
            s.quadratureOrder = obj.quadratureOrder;
            P1proc            = P1preProcessor(s);
            P1proc.preProcess();
            obj.storeParams(P1proc);
            obj.Anodal2Gauss = obj.computeA();
            lhs = obj.diffReacProb.computeLHS(obj.epsilon);
            obj.LHS = decomposition(lhs);
        end

        function x0 = getP0fromP1(obj,x)
            obj.x_reg =  obj.getP1fromP1(x);
            x0 = zeros(obj.mesh.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*obj.x_reg;
            end
        end

        function RHS = integrate_L2_function_with_shape_function(obj,x)
            RHS = obj.M*x;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                lhs = obj.diffReacProb.computeLHS(epsilon);
                obj.LHS = decomposition(lhs);
            end
        end

        function x_reg = getP1fromP1(obj,x)
            RHS = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solve_filter(RHS);
        end

        function x_reg = getP1fromP0(obj,x0)
            RHS = obj.integrateP1FunctionWithShapeFunction(x0);
            x_reg = obj.solveFilter(RHS);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.createDiffReacProblem(cParams);
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function createMassMatrix(obj)
            ss.name        = 'x';
            ss.mesh        = obj.mesh;
            s.dim          = DimensionScalar(ss);
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            lhs = LHSintegrator.create(s);
            obj.M = lhs.compute();
        end

        function storeParams(obj,P1proc)
            obj.quadrature = P1proc.quadrature;
            obj.interp     = P1proc.interp;
            obj.geometry   = P1proc.geometry;
        end

        function createDiffReacProblem(obj,cParams)
            s = cParams.femSettings;
            if isfield(cParams.femSettings,'LHStype')
                s.LHStype = cParams.femSettings.LHStype;
            else
                s.LHStype = 'DiffReactNeumann';
            end
            if isprop(cParams,'mesh')
                s.mesh = cParams.mesh;
            end
            s.type = 'DIFF-REACT';
            obj.diffReacProb = FEM.create(s);
        end

        function A_nodal_2_gauss = computeA(obj)
            s.nnode  = obj.mesh.nnode;
            s.nelem  = obj.mesh.nelem;
            s.npnod  = obj.mesh.npnod;
            s.ngaus  = obj.quadrature.ngaus;
            s.connec = obj.mesh.connec;
            s.shape  = obj.interp.shape;
            Acomp = Anodal2gausComputer(s);
            Acomp.compute();
            A_nodal_2_gauss = Acomp.A_nodal_2_gauss;
        end

        function x_reg = solve_filter(obj,RHS)
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end

        function x_reg = solveFilter(obj,RHS)
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end

        function intX = integrateP1FunctionWithShapeFunction(obj,x)
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