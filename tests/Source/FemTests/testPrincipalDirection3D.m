classdef testPrincipalDirection3D < testPrincipalDirection
    
    methods (Access = public)
       
        function obj = testPrincipalDirection3D()
            obj.init();
            obj.createTensor();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.stressDim = 6;
            obj.pdim = '3D';
            obj.nelem = 6400; 
            obj.nGaus = 3;            
        end
        
    end
    
end