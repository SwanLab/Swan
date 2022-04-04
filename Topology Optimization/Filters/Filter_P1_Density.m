classdef Filter_P1_Density < NewFilter
    
    
    methods (Access = public)
        
        function obj = Filter_P1_Density(cParams)
            obj.varInitialization(cParams);
            obj.createFilterKernel();
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
      
        function x0 = computeP0fromP1(obj,x)
            x0 = obj.Kernel*x;
        end
        
        function createFilterKernel(obj)
            P = obj.Poper.value;
            obj.Kernel = P*obj.M;
        end
        
        function intX = integrateRHS(obj,x)
            intX = zeros(obj.mesh.nelem,1);
            ngaus = size(x,2);
            for igaus = 1:ngaus
                dvolu = obj.geometry.dvolu(:,igaus);
                intX = intX + dvolu.*x(:,igaus);
            end
        end
        
    end
    
end