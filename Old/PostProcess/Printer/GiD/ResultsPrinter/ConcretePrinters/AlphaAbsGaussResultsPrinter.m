classdef AlphaAbsGaussResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'AlphaAbsGauss';
        hasGaussData = true;
    end
    
    properties (Access = private)
        fieldName = 'AlphaAbs';
    end
    
    methods (Access = public)
        
        function obj = AlphaAbsGaussResultsPrinter(d)
            obj.init(d);
        end
        
        function printResults(obj,iter,fileID)
            field(1,:,:) = obj.fields;
            dS = obj.createVectorGaussDataBase(iter,fileID,field, obj.fieldName,'OnGaussPoints','A');
            VectorGaussPrinter(dS);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.alphaAbs;
        end
        
        function createHeadPrinter(obj,d,dh)
            obj.headPrinter = GaussHeadPrinter(d,dh);
        end
        
    end
    
end