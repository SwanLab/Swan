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
            obj.createQuadrature(cParams);
            obj.createPoperator();
        end

        function xReg = compute(obj,fun,quadType) % computeWorking DELETE ASAP
            switch class(fun)
                case 'FGaussDiscontinuousFunction'
                    xReg = obj.getP1Function(fun,quadType);
                otherwise
                    xReg = obj.getFGaussFunction(fun,quadType);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
        end

        function createQuadrature(obj,cParams)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(cParams.quadType);
            obj.quadrature = q;
        end

        function createPoperator(obj)
            s.mesh    = obj.mesh;
            obj.Poper = Poperator(s);
        end

        function int = computeRHSintegrator(obj,unfFun,quadType)
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            s.mesh     = unfFun.unfittedMesh;
            int        = RHSintegrator.create(s);
        end

        function xG = expressInFilterGaussPoints(obj,x)
            x0 = zeros(length(x),obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = x;
            end
            ngaus        = obj.quadrature.ngaus;
            nelem        = obj.mesh.nelem;
            s.fValues    = reshape(x0',[1,ngaus,nelem]);
            s.mesh       = obj.mesh;
            s.quadrature = obj.quadrature;
            xG           = FGaussDiscontinuousFunction(s);
        end

        function xReg = getP1Function(obj,unfFun,quadType) % NOT OPERATIVE
            test         = P0Function.create(obj.mesh, 1);
            int          = obj.computeRHSintegrator(unfFun,quadType);
            P            = obj.Poper.value;
            A            = P';
            b            = int.integrateInDomain(unfFun,test);
            p.fValues    = A*b;
            p.mesh       = obj.mesh;
            xReg         = P1Function(p);
        end

        function xReg = getFGaussFunction(obj,unfFun,quadType)
            test         = P1Function.create(obj.mesh, 1);
            int          = obj.computeRHSintegrator(unfFun,quadType);
            P            = obj.Poper.value;
            A            = P;
            b            = int.integrateInDomain(unfFun,test);
            xR           = A*b;
            xReg         = obj.expressInFilterGaussPoints(xR);
        end

    end
    
end
