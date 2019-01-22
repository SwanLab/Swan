classdef ElasticityMicroResultsPrinter < ResultsPrinter  
    
    properties (Access = protected)
        simulationStr = 'ElasticityMicroResultsPrinter';
        hasGaussData = true;       
    end
    
    properties (Access = private)
        stressStrBase = 'Stress';
        strainStrBase = 'Strain';
        stressFlucStrBase = 'StressFluc';
        strainFlucStrBase = 'StrainFluc';
        dispStrBase = 'Displacements';
        forStrBase  = 'Forces';
        stressStr
        strainStr
        stressFlucStr
        strainFlucStr
        dispStr
        forStr

        stressCompStr = 'S';
        strainCompStr = 'E';
        displCompStr  = 'U';
        forCompStr = 'F';
        dSig
        dSigFluc
        dStr
        dStrFluc
        dDisp
        dFor
    end
    
    methods (Access = public)
        
        function obj = ElasticityMicroResultsPrinter(d)
            obj.init(d);
        end
        
        function printResults(obj,iter,fileID)
            obj.createDataBases(iter,fileID);
            VectorNodalPrinter(obj.dDisp);
            VectorNodalPrinter(obj.dFor);
            VectorGaussPrinter(obj.dSig);
            VectorGaussPrinter(obj.dStr);
            VectorGaussPrinter(obj.dSigFluc);
            VectorGaussPrinter(obj.dStrFluc);            
        end        
        
        function setStrVariablesMicroCase(obj,n)
            n = num2str(n);
            obj.stressStr     = strcat(obj.stressStrBase,n);
            obj.strainStr     = strcat(obj.strainStrBase,n);
            obj.stressFlucStr = strcat(obj.stressFlucStrBase,n);
            obj.strainFlucStr = strcat(obj.strainFlucStrBase,n);
            obj.dispStr       = strcat(obj.dispStrBase,n);
            obj.forStr        = strcat(obj.forStrBase,n);
        end
    end
        
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.fields; 
        end                        
        
    end
    
    methods (Access = private)
        
        function createDataBases(obj,iter,fileID)
            f = obj.fields;
            obj.dDisp    = obj.createVectorDataBase(iter,fileID,f.d_u,obj.dispStr,'OnNodes',obj.displCompStr);
            obj.dFor     = obj.createVectorDataBase(iter,fileID,f.fext,obj.forStr,'OnNodes',obj.forCompStr);            
            obj.dSig     = obj.createVectorGaussDataBase(iter,fileID,f.stress, obj.stressStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStr     = obj.createVectorGaussDataBase(iter,fileID,f.strain, obj.strainStr,'OnGaussPoints',obj.strainCompStr);
            obj.dSigFluc = obj.createVectorGaussDataBase(iter,fileID,f.stress_fluct, obj.stressFlucStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStrFluc = obj.createVectorGaussDataBase(iter,fileID,f.strain_fluct, obj.strainFlucStr,'OnGaussPoints',obj.strainCompStr);
        end
        
    end
    
    
end

