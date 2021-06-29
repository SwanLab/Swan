classdef Filter_P1_LevelSet <  handle %Filter_LevelSet %& Filter_P1
    
    properties (Access = private)
        mesh        
        quadratureOrder

        Poper
        
        x
        x_reg
        
        geometry
        quadrature

        domainType        
        projector
        
        interp
    end
    
    methods (Access = public)
        
        function obj = Filter_P1_LevelSet(cParams)
            obj.init(cParams);
            obj.domainType = cParams.domainType;
            obj.createQuadrature();            
            obj.createProjector();
            obj.createInterpolation();
            obj.createGeometry();
            obj.createPoperator(cParams);
            obj.disableDelaunayWarning();      
        end
        
        function preProcess(obj)
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
    
    methods (Access = protected)
        
        function x0 = computeP0fromP1(obj,x)
            xN = obj.projector.project(x);
            P  = obj.Poper.value;            
            x0 = P*xN;
        end
        
    end
    
    methods (Access = private)
        
        function createPoperator(obj,cPar)            
            cParams.nnode  = obj.mesh.nnode;
            cParams.npnod  = obj.mesh.npnod;
            cParams.connec = obj.mesh.connec;
            cParams.nelem  = obj.mesh.nelem;
            cParams.diffReactEq = cPar.femSettings;            
            obj.Poper = Poperator(cParams);
        end
        
        function createProjector(obj)
            cParams.mesh = obj.mesh;
            cParams.domainType = obj.domainType;
            cParams.type = obj.mesh.type;
            obj.projector = ShapeFunctionProjector.create(cParams);                                    
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
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end
        
        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.type);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end        
        
        function createInterpolation(obj)
            obj.interp = Interpolation.create(obj.mesh,'LINEAR');    
        end
        
        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interp);
        end

    end
    
     methods (Access = private, Static)
     
        function disableDelaunayWarning()
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
        end
        
    end
    
    
end
