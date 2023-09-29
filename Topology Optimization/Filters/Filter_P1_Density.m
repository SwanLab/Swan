classdef Filter_P1_Density < handle

    properties (Access = private)
        mesh
        Poper
        quadrature
    end

    methods (Access = public)

        function obj = Filter_P1_Density(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createPoperator();
        end

        function xReg = getP1Function(obj,f,quadType)
            s.quadType = quadType;
            s.fun      = f;
            s.trial    = P0Function.create(obj.mesh, 1);
            in         = obj.computeRHSintegrator(s);
            P          = obj.Poper.value;
            A          = P';
            b          = in.RHS;
            p.fValues  = A*b;
            p.mesh     = obj.mesh;
            xReg       = P1Function(p);
        end

        function xReg = getP0Function(obj,f,quadType)
            s.quadType = quadType;
            s.fun      = f;
            s.trial    = P1Function.create(obj.mesh, 1);
            in         = obj.computeRHSintegrator(s);
            P          = obj.Poper.value;
            A          = P;
            b          = in.RHS;
            xR         = A*b;
            x0         = zeros(length(xR),obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = xR; % Not a P0!
            end
            xReg = x0;
            % FGauss...?
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function createPoperator(obj)
            s.mesh    = obj.mesh;
            obj.Poper = Poperator(s);
        end

        function rhs = computeRHSintegrator(obj,cParams)
            s.type     = 'functionWithShapeFunction';
            s.quadType = cParams.quadType;
            s.mesh     = obj.mesh;
            s.fun      = cParams.fun;
            s.trial    = cParams.trial;
            rhs        = RHSintegrator.create(s);
        end

    end

end