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
            obj.createRHSProjector(cParams);
            obj.createPoperator();
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

        function xReg = getP1Function(obj,f,quadType)
            s.quadType = quadType;
            s.fun      = f;
            s.trial    = P0Function.create(obj.mesh, 1);
            in         = obj.computeRHSintegrator(s);
            P          = obj.Poper.value;
            A          = P';
            b          = in.RHS;
            xReg       = A*b;
        end

        function xReg = getP0Function(obj,f,quadType)
            levelSet = f.fValues;
            b        = obj.projector.project(levelSet);
            P        = obj.Poper.value;
            A        = P;
            xR       = A*b;
            x0       = zeros(length(xR),obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = xR;
            end
            xReg = x0;
        end

    end

    methods (Access = private)

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function createPoperator(obj)
            s.mesh   = obj.mesh;
            obj.Poper = Poperator(s);
        end
        
        function createRHSProjector(obj,cParams)
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
