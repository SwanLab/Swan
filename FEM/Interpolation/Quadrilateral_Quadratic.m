classdef Quadrilateral_Quadratic < Interpolation
    
    methods (Access = public)
        
        function obj = Quadrilateral_Quadratic(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.computeShapes(posgp);
            obj.computeShapeDerivatives(posgp);
        end
        
    end
    
    methods (Access = private)
        
        function computeParams(obj)
            obj.type = 'QUADRILATERAL_QUADRATIC';
            obj.ndime = 2;
            obj.nnode = 9;
            obj.pos_nodes = [0,0 ; 1 0; 0,1 ; 0.5,0 ; 0.5,0.5 ; 0,0.5];
        end
        
        function computeShapes(obj,posgp)
         ngaus = size(posgp,2);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                obj.shape(:,igaus) = [];
            end
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                obj.deriv(:,:,igaus) = [];
            end
        end

    end
    
end