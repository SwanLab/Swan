classdef  StressNormVariationInCellProblem < handle
    
    properties (SetAccess = private, GetAccess = public)
        invNormalizedStressNorm 
        invStressMax 
        pNorm      
        meshSize
    end
    
    
    properties (Access = private)
        gmsFile
        outFile
        homog      
        homogDataBase
        stressNormFun 
        stressNorm
        stressMax
    end
    
    methods (Access = public)
        
        function obj = StressNormVariationInCellProblem(d)
            obj.init(d)
            obj.createFemMatOoInputData();
            obj.createStressShapeFunction();
            obj.computeStressNorms();
            obj.computeStressMax();
            obj.computeMeshSize();
        end
        
        function printVariables(obj)
            d.stressNormFun = obj.stressNormFun;
            d.outFile       = obj.outFile;           
            printer = StressNormPrinter(d);            
            printer.print();
        end
        
    end
    
    methods (Access = private)
                         
        function init(obj,d)
            obj.gmsFile = d.gmsFile;
            obj.outFile = d.outFile; 
            nNorm = 6;
            obj.pNorm = 2.^(1:nNorm);            
        end            
        
        function createFemMatOoInputData(obj)
            oD = fullfile('Output',obj.outFile);
            fullOutFile = fullfile(oD,[obj.outFile,'.m']);
            c = GmsFile2FemMatOoFileConverter(obj.gmsFile,oD,fullOutFile);
            c.convert();
        end            
        
        function createStressShapeFunction(obj)
            d.homog   = obj.homog;  
            d.homogDB = obj.homogDataBase;
            d.outFile = obj.outFile;
            sNormCreator = StressNormShapeFuncCreator(d);
            obj.stressNormFun = sNormCreator.getStressNormShape();
        end
        
        function computeStressNorms(obj)
            nNorm = length(obj.pNorm);
            obj.stressNorm = zeros(nNorm,2);
            for ip = 1:nNorm
                obj.computeStressNorm(ip);
                obj.computeInverseNormalizedStressNorm(ip)
            end
        end
        
        function computeStressNorm(obj,ip)
            p = obj.pNorm(ip);
            sF = obj.stressNormFun;
            sF.setPnorm(p);
            sN = sF.computeCostWithFullDomain();
            obj.stressNorm(ip,1) = sN;
        end
        
        function computeInverseNormalizedStressNorm(obj,ip)
            p  = obj.pNorm(ip);        
            sN = obj.stressNorm(ip,1);
            invSn = 1/((sN)^(1/p));
            obj.invNormalizedStressNorm(ip,1) = invSn;
        end
        
        function computeStressMax(obj)
            sF = obj.stressNormFun;            
            obj.stressMax    = sF.computeMaxStressWithFullDomain();
            obj.invStressMax = 1/obj.stressMax;
        end
        
        function computeMeshSize(obj)
            sF = obj.stressNormFun;                       
            microProb  = sF.getPhysicalProblems();
            mesh       = microProb{1}.mesh;
            obj.meshSize = mesh.computeMeanCellSize();            
        end
                
    end
    
end