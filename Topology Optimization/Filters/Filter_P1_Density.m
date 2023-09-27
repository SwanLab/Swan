classdef Filter_P1_Density < Filter

    properties (Access = private)
        Poper
        x
        x_reg
        quadrature
    end

    methods (Access = public)

        function obj = Filter_P1_Density(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createPoperator(cParams);
        end

        function xReg = getP1fromP0(obj,x0)
            s.fValues = x0;
            s.mesh    = obj.mesh;
            f         = P0Function(s);
            xReg      = obj.getP1Function(f);
        end

        function xReg = getP1Function(obj,f)
            P          = obj.Poper.value;
            A          = P';
            s.type     = 'functionWithShapeFunction';
            s.quadType = 'LINEAR';
            s.mesh     = obj.mesh;
            s.fun      = f;
            s.trial    = P0Function.create(obj.mesh, 1);
            in         = RHSintegrator.create(s);
            b          = in.RHS;
            xReg       = A*b;
        end
        
        function x0 = getP0fromP1(obj,x)
            if obj.xHasChanged(x)
                s.fValues = x;
                s.mesh = obj.mesh;
                f = P1Function(s);
                xR = obj.getP0Function(f);
                x0 = zeros(length(xR),obj.quadrature.ngaus);
                for igaus = 1:obj.quadrature.ngaus
                    x0(:,igaus) = xR;
                end
            else
                x0 = obj.x_reg;
            end
            obj.updateStoredValues(x,x0);
        end

        function xReg = getP0Function(obj,f)
            P          = obj.Poper.value;
            A          = P;
            s.type     = 'functionWithShapeFunction';
            s.quadType = 'QUADRATICMASS';
            s.mesh     = obj.mesh;
            s.fun      = f;
            s.trial    = P1Function.create(obj.mesh, 1);
            in         = RHSintegrator.create(s);
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

        function itHas = xHasChanged(obj,x)
            itHas = ~isequal(x,obj.x);
        end

        function updateStoredValues(obj,x,x0)
            obj.x = x;
            obj.x_reg = x0;
        end

    end

end