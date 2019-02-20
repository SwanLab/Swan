classdef StressPdependency < handle
    
    
    properties (Access = protected)
        testName = 'AmplificatorsPdependency';
    end
    
    properties (Access = private)
        smoothDB
        nonSmoothDB
        experiments
        nExp
        expName
        outPutPath
    end
    
    methods (Access = public)
        
        function obj = StressPdependency()
            obj.init();
            obj.computeStressNormProblems();
            obj.plotStressNormInfo();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            firstPart  = fullfile( '/home','alex','Dropbox');
            secondPart = fullfile('Amplificators','Images','MicroWithHole/');
            obj.outPutPath = fullfile(firstPart,secondPart);            
            obj.experiments = {'1','2','4'};
            obj.nExp = length(obj.experiments);
        end
        
        function computeStressNormProblems(obj)
            for iexp = 1:obj.nExp
                obj.expName = obj.experiments{iexp};
                obj.computeStressNormOfSmoothProblem(iexp);
                obj.computeStressNormOfNonSmoothProblem(iexp);
            end
        end
        
        function computeStressNormOfSmoothProblem(obj,iexp)
            microCase = 'SmoothRect';
            dB = obj.smoothDB;
            dB = obj.computeStressNormProblem(microCase,iexp,dB);
            obj.smoothDB = dB;
        end
        
        function computeStressNormOfNonSmoothProblem(obj,iexp)
            microCase = 'Rect';
            dB = obj.nonSmoothDB;
            dB = obj.computeStressNormProblem(microCase,iexp,dB);
            obj.nonSmoothDB = dB;
        end
        
        function dB = computeStressNormProblem(obj,microCase,iexp,dB)
            sV = obj.createStressNormProblem(microCase);
            sV.printVariables();
            dB.invStressNorm(iexp,:) = sV.invNormalizedStressNorm;
            dB.invStressMax(iexp,:)  = sV.invStressMax;
            dB.pNorm          = sV.pNorm;
            dB.meshSize(iexp) = sV.meshSize;
        end
        
        function sV = createStressNormProblem(obj,microCase)
            d.gmsFile = obj.obtainGmsFileName(microCase);
            d.outFile = obj.obtainOutPutName(microCase);
            sV = StressNormVariationInCellProblem(d);
        end
        
        function gmsFileName = obtainGmsFileName(obj,fileName)
            path = fullfile('Input','MicroElips');
            gmsFileName = fullfile(path,[fileName,obj.expName,'.msh']);
        end
        
        function f = obtainOutPutName(obj,fileName)
            f = [obj.testName,fileName,obj.expName];
        end
        
        function plotStressNormInfo(obj)
            d.smoothDB    = obj.smoothDB;
            d.nonSmoothDB = obj.nonSmoothDB;
            d.nExp        = obj.nExp;
            d.outPath     = obj.outPutPath;
            plotter = StressNormPlotters(d);
            plotter.plot();
        end
        
    end
    
end