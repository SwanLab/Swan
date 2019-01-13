classdef ElasticityMicroResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'ElasticityMicroResultsPrinter';
    end
    
    properties (Access = private)
        stressStr = 'Stress';
        strainStr = 'Strain';
        stressFlucStr = 'StressFluc';
        strainFlucStr = 'StrainFluc';
        dispStr   = 'Displacements';
        forStr      = 'Forces';
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
        
        function obj = ElasticityMicroResultsPrinter(d,dGauss)
            obj.init(d);
            obj.storeDataBase(dGauss);                        
        end
        
        function setStrVariablesCase(obj,n)
            n = num2str(n);
            obj.stressStr     = strcat(obj.stressStr,n);
            obj.strainStr     = strcat(obj.strainStr,n);
            obj.stressFlucStr = strcat(obj.stressFlucStr,n);
            obj.strainFlucStr = strcat(obj.strainFlucStr,n);
            obj.dispStr       = strcat(obj.dispStr,n);
            obj.forStr        = strcat(obj.forStr,n);
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

