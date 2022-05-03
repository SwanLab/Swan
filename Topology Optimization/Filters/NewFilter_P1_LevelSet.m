classdef NewFilter_P1_LevelSet <  handle
    
    properties (Access = private)
        Poper
        x
        x_reg
        geometry
        quadrature
        projector
    end

    properties (Access = private)
        mesh
        quadratureOrder
    end
    
    methods (Access = public)
        
        function obj = NewFilter_P1_LevelSet(cParams)
            obj.init(cParams);
            obj.createProjector(cParams);
            obj.createPoperator(cParams);
            obj.disableDelaunayWarning();
        end
        
        function preProcess(obj,cParams)
            s.mesh = obj.mesh;
            s.quadratureOrder = obj.quadratureOrder;
            P1proc = P1preProcessor(s);
            P1proc.preProcess();
            obj.quadrature = P1proc.quadrature;
            obj.geometry = P1proc.geometry;
        end
        
        function x_reg = getP1fromP0(obj,x0)
            RHS = obj.integrateRHS(x0);
            P = obj.Poper.value;
            x_reg = P'*RHS;
        end
        
        function x0 = getP0fromP1(obj,x)
            if obj.xHasChanged(x)
                xR = obj.computeP0fromP1(x);
                x0 = zeros(length(xR),obj.quadrature.ngaus);
                for igaus = 1:obj.quadrature.ngaus
                    x0(:,igaus) = xR;
                end
            else
                x0 = obj.x_reg;
            end
            obj.updateStoredValues(x,x0);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function x0 = computeP0fromP1(obj,x)
            xN = obj.projector.project(x);
            P  = obj.Poper.value;
            x0 = P*xN;
        end

        function createPoperator(obj,cParams)
            s.nnode  = obj.mesh.nnode;
            s.npnod  = obj.mesh.npnod;
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
        
        function intX = integrateRHS(obj,x)
            intX = zeros(obj.mesh.nelem,1);
            ngaus = size(x,2);
            for igaus = 1:ngaus
                dvolu = obj.geometry.dvolu(:,igaus);
                intX = intX + dvolu.*x(:,igaus);
            end
        end
        
        function itHas = xHasChanged(obj,x)
            itHas = ~isequal(x,obj.x);
        end
        
        function updateStoredValues(obj,x,x0)
            obj.x = x;
            obj.x_reg = x0;
        end

    end
    
     methods (Access = private, Static)
     
        function disableDelaunayWarning()
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
        end
        
    end
    
    
end
