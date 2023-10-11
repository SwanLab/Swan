classdef FilterP1Unfitted <  handle
    
    properties (Access = private)
        mesh
        Poper
        quadrature
        levelSet
    end
    
    methods (Access = public)
        
        function obj = FilterP1Unfitted(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createPoperator();
        end

        function xReg = getP1Function(obj,f,quadType)
            test      = P0Function.create(obj.mesh, 1);
            int       = obj.computeRHSintegrator(quadType);
            b         = int.integrateInDomain(f,test);
            P         = obj.Poper.value;
            A         = P';
            p.fValues = A*b;
            p.mesh    = obj.mesh;
            xReg      = P1Function(p);
        end

        function xReg = getP0Function(obj,f,quadType)
            ls    = f.fValues;
            test  = P1Function.create(obj.mesh, 1);
            int   = obj.computeRHSintegrator(quadType);
            b     = int.integrateInDomain(ones(size(ls)),test); % to be included in CharFun
            P     = obj.Poper.value;
            A     = P;
            xR    = A*b;
            x0    = zeros(length(xR),obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = xR;
            end
            ngaus        = obj.quadrature.ngaus;
            nelem        = obj.mesh.nelem;
            s.fValues    = reshape(x0',[1,ngaus,nelem]);
            s.mesh       = obj.mesh;
            s.quadrature = obj.quadrature;
            xReg         = FGaussDiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.levelSet = cParams.designVariable;
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

        function int = computeRHSintegrator(obj,quadType)
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            s.mesh     = obj.levelSet.getUnfittedMesh();
            int        = RHSintegrator.create(s);
        end

    end
    
end
