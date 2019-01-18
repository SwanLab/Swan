classdef ElasticityResultsPrinter < ResultsPrinter ...
        & GaussResultsPrinter
    
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
        
        function obj = ElasticityResultsPrinter(d)
            obj.init(d);
        end
        
        function setStrVariablesNames(obj,s1,s2,s3)
            obj.stressStr = s1;
            obj.strainStr = s2;
            obj.dispStr   = s3;
        end
        
        function storeResultsInfo(obj,d)
            obj.storeQuadInfo(d);
            obj.fields = d.variables;
        end
        
        
    end
    
    methods (Access = protected)
        
        
        function printHeader(obj)
            d.fileID = obj.fileID;
            d.etype = obj.etype;
            d.ndim  = obj.ndim;
            d.gaussDescriptor = obj.gaussDescriptor;
            d.posgp = obj.posgp;
            d.ngaus = obj.ngaus;
            obj.headPrinter.print(d);
        end
        
        function printResults(obj)
            obj.createDataBases(obj.fields);
            VectorNodalPrinter(obj.dV);
            VectorGaussPrinter(obj.dSig);
            VectorGaussPrinter(obj.dStr);
        end
        
    end
    
    methods (Access = private)
        
        function createDataBases(obj,fields)
            f = fields;
            obj.dV   = obj.createVectorDataBase(f.d_u,obj.dispStr,'OnNodes',obj.displCompStr);
            obj.dSig = obj.createVectorGaussDataBase(f.stress, obj.stressStr,'OnGaussPoints',obj.stressCompStr);
            obj.dStr = obj.createVectorGaussDataBase(f.strain, obj.strainStr,'OnGaussPoints',obj.strainCompStr);
        end
        
    end
    
    
end

