classdef FunctionResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'FeFunction';
        hasGaussData = true;
    end
    
    properties (Access = private)
        valuesStr = 'FValues';
        strainStr = 'Strain';
        dispStr   = 'Displacements';
        stressCompStr = 'S';
        strainCompStr = 'E';
        displCompStr  = 'U';
        fvalsCompStr  = 'V';
        dSig
        dStr
        dV
        nameCase
        fVals
    end
    
    methods (Access = public)
        
        function obj = FunctionResultsPrinter(d)
            obj.init(d);
            obj.nameCase = d.name;
            obj.setStrVariablesNames()
        end
        
        function printResults(obj,iter,fileID)
            obj.createDataBases(obj.fields,iter,fileID);
            VectorGaussPrinter(obj.fVals);
%             VectorGaussPrinter(obj.dStr);
%             VectorNodalPrinter(obj.dV);
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
%             obj.stressStr = [obj.stressStr,n];
%             obj.strainStr = [obj.strainStr,n];
%             obj.dispStr   = [obj.dispStr,n];
        end
        
        function createDataBases(obj,fields,iter,fileID)
            f = fields;
            obj.fVals = obj.createVectorGaussDataBase(iter,fileID,f{1}, obj.valuesStr,'OnGaussPoints',obj.fvalsCompStr);
%             obj.dV   = obj.createVectorDataBase(iter,fileID,f.u,obj.dispStr,'OnNodes',obj.displCompStr);
%             obj.dSig = obj.createVectorGaussDataBase(iter,fileID,f.stress, obj.stressStr,'OnGaussPoints',obj.stressCompStr);
%             obj.dStr = obj.createVectorGaussDataBase(iter,fileID,f.strain, obj.strainStr,'OnGaussPoints',obj.strainCompStr);
        end
        
    end
    
end

