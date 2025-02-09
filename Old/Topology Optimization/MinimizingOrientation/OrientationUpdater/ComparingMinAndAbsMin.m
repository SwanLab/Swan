classdef ComparingMinAndAbsMin < handle
    
    properties (Access = private)
       orientationUpdaterType
       orientationProblem
       nStressFigure
       cost
       icase
       cases
       resultsDir
    end
    
    methods (Access = public)
        
        function obj = ComparingMinAndAbsMin()
            obj.resultsDir = '/home/alex/Dropbox/Amplificators/GregoireMeeting7/';
            obj.cases = {'MinimumEigenValue',...
                          'MaximumEigenValue',...
                          'MinimumAbsEigenValue',...
                          'MaximumAbsEigenValue'};
            obj.computeOrientations();
            obj.plotAllCosts();
        end
        
    end
    
    methods (Access = private)
        
        function computeOrientations(obj)
            for ic = 1:numel(obj.cases)
                obj.icase = ic;
                obj.computeOrientationCase();
            end
        end
        
        function computeOrientationCase(obj)
            obj.nStressFigure = obj.icase;
            obj.orientationUpdaterType = obj.cases{obj.icase};
            obj.computeOrientationProblem();
            obj.plotStress();
            obj.cost{obj.icase} = obj.orientationProblem.cost;
        end
        
        function computeOrientationProblem(obj)
            cParams.plotting = true;
            cParams.orientationUpdaterType = obj.orientationUpdaterType;
            obj.orientationProblem = OrientationAsDesignVariable(cParams);
            obj.orientationProblem.compute();
        end
        
        function plotStress(obj)
            obj.orientationProblem.nStressFigure = obj.nStressFigure;
            fID = obj.orientationProblem.plotStress();
            figName = [obj.resultsDir,obj.cases{obj.icase}];
            title(obj.cases{obj.icase})
            print(fID,figName,'-dpng');
        end

        function plotAllCosts(obj)
           f = figure;
           hold on;
           for ic = 1:numel(obj.cases)
               c = obj.cost{ic};
               h{ic} = plot(c,'-+');
           end
           legend(obj.cases)
           pp = plotPrinter(f,h);
           pp.print([obj.resultsDir,'CostComparison'])
        end
        
    end
    
end

