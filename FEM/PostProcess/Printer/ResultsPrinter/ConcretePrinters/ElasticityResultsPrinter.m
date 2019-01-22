classdef ElasticityResultsPrinter < ResultsPrinter ...
                                    
    
    
    properties (Access = protected)
        simulationStr = 'ElasticityResults';
        hasGaussData = true;       
    end
    
    properties (Access = private)
        stressStr = 'Stress';
        strainStr = 'Strain';
        dispStr   = 'Displacements';
        stressCompStr = 'S';
        strainCompStr = 'E';
        displCompStr  = 'U';
        dSig
        dStr
        dV
    end
    
    methods (Access = public)
        
        function obj = ElasticityResultsPrinter(d)
            obj.init(d);
        end
        
        function printResults(obj,iter,fileID)
            obj.createDataBases(obj.fields,iter,fileID);
            VectorNodalPrinter(obj.dV);
            VectorGaussPrinter(obj.dSig);
            VectorGaussPrinter(obj.dStr);
        end
        
        function setStrVariablesNames(obj,s1,s2,s3)
            obj.stressStr = s1;
            obj.strainStr = s2;
            obj.dispStr   = s3;
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.fields; 
        end
                    
    end
    
    methods (Access = private)
        
        function createDataBases(obj,fields,iter,fileID)
            f = fields;
            obj.dV   = obj.createVectorDataBase(iter,fileID,f.d_u,obj.dispStr,'OnNodes',obj.displCompStr);
            obj.dSig = obj.createVectorGaussDataBase(iter,fileID,f.stress, obj.stressStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStr = obj.createVectorGaussDataBase(iter,fileID,f.strain, obj.strainStr,'OnGaussPoints',obj.strainCompStr);
        end
        
    end
    
    
end

