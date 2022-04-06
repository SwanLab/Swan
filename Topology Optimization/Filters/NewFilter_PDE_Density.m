classdef NewFilter_PDE_Density < handle
    
    properties (Access = private)
        Anodal2Gauss
        LHS
        epsilon
        mesh
        quadratureOrder
        diffReacProb
        nelem
        ngaus
        npnod
        nnode
        shape
        quadrature
        geometry
        interp
    end

    methods (Access = public)
        
        function obj = NewFilter_PDE_Density(cParams)
            obj.init(cParams);
            obj.epsilon = cParams.mesh.computeMeanCellSize();
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
            RHS = obj.diffReacProb.getM()*x;
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

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.createDiffReacProblem(cParams);
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function A_nodal_2_gauss = computeA(obj)
            A0 = sparse(obj.nelem,obj.npnod);
            A_nodal_2_gauss = cell(obj.ngaus,1);
            fn = ones(1,obj.npnod);

            nodes = obj.mesh.connec;
            fg = zeros(obj.ngaus,obj.nelem);

            for igaus = 1:obj.ngaus
                A_nodal_2_gauss{igaus} = A0;
                for inode = 1:obj.nnode
                    node   = nodes(:,inode);
                    fN     = fn(node);
                    shapeN = obj.shape(inode,igaus);
                    fg(igaus,:) = fg(igaus,:) + shapeN*fN;
                    Ni = ones(obj.nelem,1)*shapeN;
                    A  = sparse(1:obj.nelem,node,Ni,obj.nelem,obj.npnod);
                    A_nodal_2_gauss{igaus} = A_nodal_2_gauss{igaus} + A;
                end
            end
        end

        function computeLHS(obj)
            lhs = obj.diffReacProb.computeLHS();
            obj.LHS = decomposition(lhs);
        end

        function createDiffReacProblem(obj,cParams)
            s = cParams.femSettings;
            s.mesh = cParams.mesh;
            switch s.scale
                case 'MACRO'
                    obj.diffReacProb = NewDiffReactProblem(s);
                case 'MICRO'
                    obj.diffReacProb = NewDiffReactProblemMicro(s);
            end
        end

        function storeParams(obj,P1proc)
            obj.quadrature = P1proc.quadrature;
            obj.interp     = P1proc.interp;
            obj.geometry   = P1proc.geometry;
            obj.nelem      = obj.mesh.nelem;
            obj.nnode      = obj.mesh.nnode;
            obj.npnod      = obj.mesh.npnod;
            obj.ngaus      = obj.quadrature.ngaus;
            obj.shape      = obj.interp.shape;
        end

        function x_reg = solve_filter(obj,RHS)
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end

    end

end