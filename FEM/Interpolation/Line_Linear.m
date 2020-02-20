classdef Line_Linear < Interpolation
    
    methods (Access = public)
        
        function obj = Line_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
            obj.computeCases();            
            obj.computeCoordAndConnec();
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.computeShapes(posgp);
            obj.computeShapeDerivatives(posgp);
        end
        
    end
    
    methods (Access = private)

        function computeParams(obj)
            obj.type = obj.mesh.geometryType;
            obj.ndime = 1;
            obj.nnode = 2;
            obj.pos_nodes = [-1; 1];
            obj.isoDv = 2;  
        end
        
        function computeShapes(obj,posgp)
            ngaus = size(posgp,2);
            s = posgp(1,:);
            I = ones(1,ngaus);
            obj.shape = [I-s;s+1]/2;            
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);            
            obj.deriv = repmat([-0.5,0.5],1,1,ngaus);
        end
        
        function computeCases(obj)
            obj.iteration = [1; 
                             2 ];                            
        end
        
    end
    
end
