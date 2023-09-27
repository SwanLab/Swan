classdef Filter_P1_Density < Filter

    properties (Access = private)
        Poper
        quadrature
    end

    methods (Access = public)

        function obj = Filter_P1_Density(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createPoperator(cParams);
        end

        function rhs = computeRHSintegrator(obj,cParams)
            s.type     = 'functionWithShapeFunction';
            s.quadType = cParams.quadType;
            s.mesh     = obj.mesh;
            s.fun      = cParams.fun;
            s.trial    = cParams.trial;
            rhs        = RHSintegrator.create(s);
        end

        function xReg = getP1fromP0(obj,x0)
            nelem     = size(x0,1);
            ngaus     = size(x0,2);
            s.fValues = reshape(x0',[1,ngaus,nelem]);
            s.mesh    = obj.mesh;
            s.quadrature = obj.quadrature;
            f         = FGaussDiscontinuousFunction(s);
            xReg      = obj.getP1Function(f);
        end

        function xReg = getP1Function(obj,f)
            s.quadType = 'LINEAR';
            s.fun      = f;
            s.trial    = P0Function.create(obj.mesh, 1);
            in         = obj.computeRHSintegrator(s);
            P          = obj.Poper.value;
            A          = P';
            b          = in.RHS;
            xReg       = A*b;
        end

        function x0 = getP0fromP1(obj,x)
            s.fValues = x;
            s.mesh = obj.mesh;
            f = P1Function(s);
            xR = obj.getP0Function(f);
            x0 = zeros(length(xR),obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = xR;
            end
        end

        function xReg = getP0Function(obj,f)
            s.quadType = 'QUADRATICMASS';
            s.fun      = f;
            s.trial    = P1Function.create(obj.mesh, 1);
            in         = obj.computeRHSintegrator(s);
            P          = obj.Poper.value;
            A          = P;
            b          = in.RHS;
            xReg       = A*b;
        end

    end

    methods (Access = private)

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function createPoperator(obj,cPar)
            cParams.nelem  = obj.mesh.nelem;
            cParams.nnode  = obj.mesh.nnodeElem;
            cParams.npnod  = obj.mesh.nnodes;
            cParams.connec = obj.mesh.connec;
            cParams.diffReactEq = cPar.femSettings;
            obj.Poper = Poperator(cParams);
        end

    end

end