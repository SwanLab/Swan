classdef Filter_P1_LevelSet <  Filter
    
    properties (Access = private)
        Poper
        projector
        quadrature
    end
    
    methods (Access = public)
        
        function obj = Filter_P1_LevelSet(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createProjector(cParams);
            obj.createPoperator(cParams);
            obj.disableDelaunayWarning();
        end

        function rhs = computeRHSintegrator(obj,cParams)
            s.type     = 'functionWithShapeFunction';
            s.quadType = cParams.quadType;
            s.mesh     = obj.mesh;
            s.fun      = cParams.fun;
            s.trial    = cParams.trial;
            rhs        = RHSintegrator.create(s);
        end

        function x_reg = getP1fromP0(obj,x0)
            nelem     = size(x0,1);
            ngaus     = size(x0,2);
            s.fValues = reshape(x0',[1,ngaus,nelem]);
            s.mesh = obj.mesh;
            s.quadrature = obj.quadrature;
            f = FGaussDiscontinuousFunction(s);
            x_reg = obj.getP1Function(f);
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
                xR = obj.computeP0fromP1(x);
                x0 = zeros(length(xR),obj.quadrature.ngaus);
                for igaus = 1:obj.quadrature.ngaus
                    x0(:,igaus) = xR;
                end
        end

        % getP0Function ... time to include the unfitted mesh!

    end

    methods (Access = private)

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function x0 = computeP0fromP1(obj,x)
            xN = obj.projector.project(x);
            P  = obj.Poper.value;
            x0 = P*xN;
        end

        function createPoperator(obj,cParams)
            s.nnode  = obj.mesh.nnodeElem;
            s.npnod  = obj.mesh.nnodes;
            s.connec = obj.mesh.connec;
            s.nelem  = obj.mesh.nelem;
            s.diffReactEq = cParams.femSettings;
            obj.Poper = Poperator(s);
        end
        
        function createProjector(obj,cParams)
            s.mesh = cParams.mesh;
            s.type = cParams.mesh.type;
            s.domainType = cParams.domainType;
            obj.projector = ShapeFunctionProjector.create(s);
        end

    end
    
     methods (Access = private, Static)
     
        function disableDelaunayWarning()
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
        end
        
    end
    
    
end
