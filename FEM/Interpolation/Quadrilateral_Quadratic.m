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
            obj.pos_nodes = [0,0 ; 1 0 ; 0,1 ; 1,1 ; 0.5,0 ; 1,0.5 ; 0.5,1 ; 0,0.5 ; 0.5,0.5];
        end
        
        function computeShapes(obj,posgp)
         ngaus = size(posgp,2);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                obj.shape(:,igaus) = [s.*-3-t.*3+s.^2.*t.^2.*4+s.*t.*9-s.*t.^2.*6-s.^2.*t.*6+s.^2.*2+t.^2.*2+1;
                    -s + s.^2.*t.^2.*4 + s.*t.*3 - s.*t.^2.*2 - s.^2.*t.*6 + s.^2.*2;
                    -t + s.^2.*t.^2.*4 + s.*t.*3 - s.*t.^2.*6 - s.^2.*t.*2 + t.^2.*2;
                    s.^2.*t.^2.*4 + s.*t - s.*t.^2.*2 - s.^2.*t.*2;
                    s.*4 - s.^2.*t.^2.*8 - s.*t.*12 + s.*t.^2.*8 + s.^2.*t.*12 - s.^2.*4;
                    s.^2.*t.^2.*-8 - s.*t.*4 + s.*t.^2.*4 + s.^2.*t.*8;
                    s.^2.*t.^2.*-8 - s.*t.*4 + s.*t.^2.*8 + s.^2.*t.*4;
                    t.*4 - s.^2.*t.^2.*8 - s.*t.*12 + s.*t.^2.*12 + s.^2.*t.*8 - t.^2.*4;
                    s.^2.*t.^2.*16 + 1 + s.*t.*16 - s.*t.^2.*16 - s.^2.*t.*16];
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