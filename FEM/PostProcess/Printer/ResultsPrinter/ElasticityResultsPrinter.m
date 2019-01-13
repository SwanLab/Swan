classdef ElasticityResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'ElasticityResults';
    end
    
    properties (Access = private)
        stressStr = 'Stress';
        strainStr = 'Strain';
        dispStr   = 'Displacements';
        stressCompStr = 'S';
        strainCompStr = 'E';
        displCompStr  = 'U';
        headPrinter = GaussHeadPrinter;
        dSig
        dStr
        dV
    end
    
    methods (Access = public)
        
        function obj = ElasticityResultsPrinter(d,dGauss)
            obj.init(d);
            obj.storeDataBase(dGauss);            
        end
        
        function setStrVariablesNames(obj,s1,s2,s3)
            obj.stressStr = s1;
            obj.strainStr = s2;
            obj.dispStr   = s3;
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
            VectorNodalPrinter(obj.dV);
            VectorGaussPrinter(obj.dSig);
            VectorGaussPrinter(obj.dStr);
        end
        
    end
    
    methods (Access = private)
        
        function createDataBases(obj)
            f = obj.fields;
            obj.dV   = obj.createVectorDataBase(f.d_u,obj.dispStr,'OnNodes',obj.displCompStr);
            obj.dSig = obj.createVectorGaussDataBase(f.stress, obj.stressStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStr = obj.createVectorGaussDataBase(f.strain, obj.strainStr,'OnGaussPoints',obj.strainCompStr);
        end
        
    end
    
    
end

