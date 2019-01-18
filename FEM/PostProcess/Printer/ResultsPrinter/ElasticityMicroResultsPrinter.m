classdef ElasticityMicroResultsPrinter < ResultsPrinter ...
        & GaussResultsPrinter 
    
    properties (Access = protected)
        simulationStr = 'ElasticityMicroResultsPrinter';
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
        headPrinter = GaussHeadPrinter;
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
        
        function setStrVariablesMicroCase(obj,n)
            n = num2str(n);
            obj.stressStr     = strcat(obj.stressStrBase,n);
            obj.strainStr     = strcat(obj.strainStrBase,n);
            obj.stressFlucStr = strcat(obj.stressFlucStrBase,n);
            obj.strainFlucStr = strcat(obj.strainFlucStrBase,n);
            obj.dispStr       = strcat(obj.dispStrBase,n);
            obj.forStr        = strcat(obj.forStrBase,n);
        end
        
        function storeResultsInfo(obj,d)
            obj.storeQuadInfo(d);
            obj.fields = d.variables;
        end
        
    end
    
    methods (Access = protected)
        
        function printHeader(obj)
            d.fileID = obj.fileID;
            d.gaussDescriptor = obj.gaussDescriptor;
            d.etype = obj.etype;
            d.ngaus = obj.ngaus;
            d.ndim  = obj.ndim;
            d.posgp = obj.posgp;
            obj.headPrinter.print(d);
        end
        
        function printResults(obj)
            obj.createDataBases();
            VectorNodalPrinter(obj.dDisp);
            VectorNodalPrinter(obj.dFor);
            VectorGaussPrinter(obj.dSig);
            VectorGaussPrinter(obj.dStr);
            VectorGaussPrinter(obj.dSigFluc);
            VectorGaussPrinter(obj.dStrFluc);            
        end
        
    end
    
    methods (Access = private)
        
        function createDataBases(obj)
            f = obj.fields;
            obj.dDisp    = obj.createVectorDataBase(f.d_u,obj.dispStr,'OnNodes',obj.displCompStr);
            obj.dFor     = obj.createVectorDataBase(f.fext,obj.forStr,'OnNodes',obj.forCompStr);            
            obj.dSig     = obj.createVectorGaussDataBase(f.stress, obj.stressStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStr     = obj.createVectorGaussDataBase(f.strain, obj.strainStr,'OnGaussPoints',obj.strainCompStr);
            obj.dSigFluc = obj.createVectorGaussDataBase(f.stress_fluct, obj.stressFlucStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStrFluc = obj.createVectorGaussDataBase(f.strain_fluct, obj.strainFlucStr,'OnGaussPoints',obj.strainCompStr);
        end
        
    end
    
    
end

