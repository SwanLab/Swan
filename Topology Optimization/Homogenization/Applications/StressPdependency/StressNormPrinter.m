classdef StressNormPrinter < handle
    
    properties (Access = private)
        postProcess
        postProcessDB
        stressNormFun
        outFile
        iter = 0;
    end
    
    methods (Access = public)
        
        function obj = StressNormPrinter(d)
            obj.init(d);
            obj.createPostProcessDataBase();
            obj.createPostProcess();
        end
        
        function print(obj)
            sF = obj.stressNormFun;
            microProb  = sF.getPhysicalProblems();
            d.quad     = microProb{1}.element.quadrature;
            d.fields   = microProb{1}.variables;
            obj.postProcess.print(obj.iter,d);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.stressNormFun = d.stressNormFun;
            obj.outFile       = d.outFile;
        end
        
        function createPostProcess(obj)
            postCase = 'ElasticityMicro';
            obj.postProcess = Postprocess(postCase,obj.postProcessDB);
        end
        
        function createPostProcessDataBase(obj)
            sF         = obj.stressNormFun;
            microProb  = sF.getPhysicalProblems();
            dB.mesh    = microProb{1}.mesh;
            dB.outName = obj.outFile;
            ps = PostProcessDataBaseCreator(dB);
            obj.postProcessDB = ps.getValue();
        end
        
    end
    
end