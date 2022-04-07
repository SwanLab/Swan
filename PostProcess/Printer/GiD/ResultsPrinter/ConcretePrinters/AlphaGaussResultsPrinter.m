classdef AlphaGaussResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'AlphaGauss';
        hasGaussData = true;
    end
    
    properties (Access = private)
        fieldName = 'Alpha';
    end
    
    methods (Access = public)
        
        function obj = AlphaGaussResultsPrinter(d)
            obj.init(d);
        end
        
        function printResults(obj,iter,fileID)
            alpha(1,:,:) = obj.fields;
            dS = obj.createVectorGaussDataBase(iter,fileID,alpha, obj.fieldName,'OnGaussPoints','A');
            VectorGaussPrinter(dS);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.alpha;
        end
        
        function createHeadPrinter(obj,d,dh)
            obj.headPrinter = GaussHeadPrinter(d,dh);
        end
        
    end
    
end