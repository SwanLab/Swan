classdef Line_Linear < Interpolation
    
    methods (Access = public)
        
        function obj = Line_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
            obj.computeCases();            
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.computeShapes(posgp);
            obj.computeShapeDerivatives(posgp);
        end
        
    end
    
    methods (Access = private)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 2;
            obj.pos_nodes = [-1; 1];
            obj.isoDv = 2;  
        end
        
        function computeShapes(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);            
            s = posgp(1,:,:);
            I = ones(size(s));                 
            obj.shape = zeros(obj.nnode,ngaus,nelem);
            obj.shape(1,:,:) = 0.5*(I-s);
            obj.shape(2,:,:) = 0.5*(s+I);                 
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);                          
            obj.deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            obj.deriv(1,1,:,:) = -0.5;
            obj.deriv(1,2,:,:) = 0.5;
        end
        
        function computeCases(obj)
            obj.iteration = [1; 
                             2 ];                            
        end
        
    end
    
end
