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
        nameCase
    end
    
    methods (Access = public)
        
        function obj = ElasticityResultsPrinter(d)
            obj.init(d);
            obj.nameCase = d.name;
            obj.setStrVariablesNames()
        end
        
        function printResults(obj,iter,fileID)
            obj.createDataBases(obj.fields,iter,fileID);
            VectorGaussPrinter(obj.dSig);
            VectorGaussPrinter(obj.dStr);
            VectorNodalPrinter(obj.dV);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.fields;
        end
        
    end
    
    methods (Access = private)
        
        function setStrVariablesNames(obj)
            n = obj.nameCase;
            obj.stressStr = [obj.stressStr,n];
            obj.strainStr = [obj.strainStr,n];
            obj.dispStr   = [obj.dispStr,n];
        end
        
        function createDataBases(obj,fields,iter,fileID)
            f = fields;
            obj.dV   = obj.createVectorDataBase(iter,fileID,f.u,obj.dispStr,'OnNodes',obj.displCompStr);
            obj.dSig = obj.createVectorGaussDataBase(iter,fileID,f.stress, obj.stressStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStr = obj.createVectorGaussDataBase(iter,fileID,f.strain, obj.strainStr,'OnGaussPoints',obj.strainCompStr);
        end
        
    end
    
end

