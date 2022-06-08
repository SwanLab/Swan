classdef Filter_P1_Density < handle
    
    properties (Access = private)
        mesh
        quadratureOrder
        
        Poper

        x
        x_reg

        geometry
        quadrature
        
        M
        Kernel
        interp

        field
    end
    
    properties
        
    end
    
    methods (Access = public)
        
        function obj = Filter_P1_Density(cParams)
            obj.init(cParams);
            obj.createField();
            obj.createMassMatrix();
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.createPoperator(cParams);
            obj.createFilterKernel();
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
    
    methods (Access = private)
        
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
      
        function x0 = computeP0fromP1(obj,x)
            x0 = obj.Kernel*x;
%             s.connec = obj.mesh.connec;
%             s.type = obj.mesh.type;
%             s.fNodes = x;
%             fe = FeFunction(s);
%             x0 = fe.computeValueInCenterElement();
        end
        
        function createFilterKernel(obj)
            P = obj.Poper.value;
            obj.Kernel = P*obj.M;
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.field = Field(s);
        end
        
        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            LHS = LHSintegrator.create(s);
            %M = LHS.compute();
            obj.M = LHS.compute();
        end
        
        function intX = integrateRHS(obj,x)
            intX = zeros(obj.mesh.nelem,1);
            ngaus = size(x,2);
            for igaus = 1:ngaus
                dvolu = obj.geometry.dvolu(:,igaus);
                intX = intX + dvolu.*x(:,igaus);
            end
        end
        
        function createPoperator(obj,cPar)
            cParams.nelem  = obj.mesh.nelem;
            cParams.nnode  = obj.mesh.nnodeElem;
            cParams.npnod  = obj.mesh.nnodes;
            cParams.connec = obj.mesh.connec;
            cParams.diffReactEq = cPar.femSettings;
            obj.Poper = Poperator(cParams);
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
    
end
