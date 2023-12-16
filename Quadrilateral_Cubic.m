classdef Quadrilateral_Cubic < Interpolation
    
    methods (Access = public)
        
        function obj = Quadrilateral_Cubic(cParams)
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
            obj.pos_nodes = [0,0 ; 1 0 ; 0,1 ; 1,1 ; 0.5,0 ; 1,0.5 ; 0.5,1 ; 0,0.5 ; 0.5,0.5];
        end
        
        function computeShapes(obj,posgp)
         ngaus = size(posgp,2);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                obj.shape(:,igaus) = [4*s^2*t^2 - 6*s^2*t + 2*s^2 - 6*s*t^2 + 9*s*t - 3*s + 2*t^2 - 3*t + 1;
                    4*s^2*t^2 - 6*s^2*t + 2*s^2 - 2*s*t^2 + 3*s*t - s;
                    4*s^2*t^2 - 2*s^2*t - 6*s*t^2 + 3*s*t + 2*t^2 - t;
                    4*s^2*t^2 - 2*s^2*t - 2*s*t^2 + s*t;
                    - 8*s^2*t^2 + 12*s^2*t - 4*s^2 + 8*s*t^2 - 12*s*t + 4*s;
                    - 8*s^2*t^2 + 8*s^2*t + 4*s*t^2 - 4*s*t;
                    - 8*s^2*t^2 + 4*s^2*t + 8*s*t^2 - 4*s*t;
                    - 8*s^2*t^2 + 8*s^2*t + 12*s*t^2 - 12*s*t - 4*t^2 + 4*t;
                    16*s^2*t^2 - 16*s^2*t - 16*s*t^2 + 16*s*t];
            end
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                obj.deriv(:,:,igaus) = [
                    4*s + 9*t - 12*s*t + 8*s*t^2 - 6*t^2 - 3,
                    4*s + 3*t - 12*s*t + 8*s*t^2 - 2*t^2 - 1,
                    3*t - 4*s*t + 8*s*t^2 - 6*t^2,
                    t - 4*s*t + 8*s*t^2 - 2*t^2,
                    24*s*t - 12*t - 8*s - 16*s*t^2 + 8*t^2 + 4,
                    16*s*t - 4*t - 16*s*t^2 + 4*t^2,
                    8*s*t - 4*t - 16*s*t^2 + 8*t^2,
                    16*s*t - 12*t - 16*s*t^2 + 12*t^2,
                    16*t - 32*s*t + 32*s*t^2 - 16*t^2;
                    
                    9*s + 4*t - 12*s*t + 8*s^2*t - 6*s^2 - 3,
                    3*s - 4*s*t + 8*s^2*t - 6*s^2,
                    3*s + 4*t - 12*s*t + 8*s^2*t - 2*s^2 - 1,
                    s - 4*s*t + 8*s^2*t - 2*s^2,
                    16*s*t - 12*s - 16*s^2*t + 12*s^2,
                    8*s*t - 4*s - 16*s^2*t + 8*s^2,
                    16*s*t - 4*s - 16*s^2*t + 4*s^2,
                    24*s*t - 8*t - 12*s - 16*s^2*t + 8*s^2 + 4,
                    16*s - 32*s*t + 32*s^2*t - 16*s^2];
            end
        end

    end
    
end