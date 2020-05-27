classdef testPrincipalDirection2D < testPrincipalDirection
    
    methods (Access = public)
       
        function obj = testPrincipalDirection2D()
            obj.init();
            obj.createTensor();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.stressDim = 3;
            obj.pdim = '2D';
            obj.nelem = 6400; 
            obj.nGaus = 3;            
        end
        
    end
    
end