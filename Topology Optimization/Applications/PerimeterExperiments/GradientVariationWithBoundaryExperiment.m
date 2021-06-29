classdef GradientVariationWithBoundaryExperiment < handle
    
    properties (Access = private)
        errorGradient
        meanGradient
        devGradient
        legendPlot  
        epsilonMeshSizeRatio        
        iMesh
        gExperiment
        circleMesh
        curvature
        kappaLegend
    end
    
    properties (Access = private)
        nameCase
        inputFiles
        outputFolder
        circleCase
        levelSetParams
    end
    
    methods (Access = public)
        
        function obj = GradientVariationWithBoundaryExperiment(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            for im = 1:numel(obj.inputFiles)
                obj.iMesh = im;
                obj.createGradientVariationExperiment();
                obj.createCircleMesh();
                obj.computeGradientVariationInTheBoundary();
                obj.computeEpsilonMeshSizeRatio();
                obj.computeLegend();                
            end       
            obj.plotDeviation();
            obj.plotMean();
            obj.plotError();
        end   
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nameCase     = cParams.nameCase;
            obj.inputFiles   = cParams.inputFiles; 
            obj.outputFolder = cParams.outputFolder;
            obj.curvature    = cParams.curvature;
            obj.circleCase   = cParams.circleCase;
            obj.levelSetParams = cParams.levelSetParams;
        end
    
        function computeGradientVariationInTheBoundary(obj)
            filePlotName = [obj.outputFolder,obj.inputFiles{obj.iMesh}];
            s.filePlotName         = filePlotName;
            s.backgroundMesh       = obj.gExperiment.backgroundMesh;
            s.levelSet             = obj.gExperiment.levelSet;
            s.regularizedPerimeter = obj.gExperiment.regularizedPerimeter;
            s.domainLength         = obj.gExperiment.domainLength;            
            s.circleMesh           = obj.circleMesh;
            s.circleCase           = obj.circleCase;
            s.curvature            = obj.curvature;
            g = GradientVariationWithBoundaryComputer(s);
            g.compute();
            obj.obtainVariables(g);
        end
        
        function createCircleMesh(obj)
            s.backgroundMesh = obj.gExperiment.backgroundMesh;
            s.boundaryMesh   = obj.gExperiment.boundaryMesh;
            uMesh = UnfittedMesh(s);
            uMesh.compute(obj.gExperiment.levelSet);
            switch obj.circleCase
                case 'interior'
                    obj.circleMesh = uMesh.boundaryCutMesh;   
                    obj.kappaLegend = '$1/R$';
                case 'exterior'
                    obj.circleMesh = obj.gExperiment.boundaryMesh{1};                    
                    obj.kappaLegend = '$1/R$';                     
            end
        end
       
        function createGradientVariationExperiment(obj)
            s.inputFile = obj.inputFiles{obj.iMesh};
            s.iMesh     = obj.iMesh;
            s.levelSetParams = obj.levelSetParams;
            g = GradientVariationExperiment(s);
            obj.gExperiment = g;
        end
        
       function obtainVariables(obj,gComputer)
           g = gComputer;
           obj.errorGradient(obj.iMesh,:) = g.error;
           obj.meanGradient(obj.iMesh,:) = g.mean;
           obj.devGradient(obj.iMesh,:) = g.desv;
       end
       
        function computeLegend(obj)
            h = obj.gExperiment.backgroundMesh.computeMeanCellSize;
            L = obj.gExperiment.domainLength;
            [n,d] = rat(h/L);
            if n == 1
                hStr = [num2str(n),'/',num2str(d)];
            else
                hStr = num2str(round(h/L,3));
            end            
            obj.legendPlot{obj.iMesh} = ['$h = ',hStr,'$'];
        end
                
        function computeEpsilonMeshSizeRatio(obj)
            h = obj.gExperiment.backgroundMesh.computeMeanCellSize;
            e = obj.gExperiment.regularizedPerimeter.epsilons;
            eh = e/h;
            obj.epsilonMeshSizeRatio(obj.iMesh,:) = eh;  
        end       
       
        function plotDeviation(obj)
            name = 'DesviationCurvatureCircleInclusion';
            leg  = obj.legendPlot;
            var  = (obj.devGradient);
            [f,hP] = obj.plotVariable(leg,var);  
            obj.printPlot(f,hP,name);            
        end
        
        function plotMean(obj)
            name = 'MeanCurvatureCircleInclusion';
            leg = obj.legendPlot;
            leg{end+1} = obj.kappaLegend;
            var = obj.meanGradient;
            var(end+1,:) = obj.curvature*ones(size(var(1,:)));            
            [f,hP] = obj.plotMeanVariable(leg,var);  
            obj.printPlot(f,hP,name);            
        end    
        
        function plotError(obj)
            name = 'ErrorCurvatureCircleInclusion';
            leg  = obj.legendPlot;
            var  = (obj.errorGradient);
            [f,hP] = obj.plotVariable(leg,var);  
            obj.printPlot(f,hP,name);            
        end               
       
       function [f,hP] = plotVariable(obj,legendPlot,var)            
            x = log2(obj.epsilonMeshSizeRatio(1,:));
            y = var';
            f = figure();
            hP{1} = semilogy(x,y,'+-');
            set(gca, 'XTickLabel',[])
            set(gca,'XTick',x,'XTickLabel',{'1' '2' '4' '8' '16' '32'})
            set(gca, 'Xdir', 'reverse');
            legObj = legend(legendPlot);            
            set(legObj,'Interpreter','latex','Location','best');
            xlabel('$\varepsilon/h$','interpreter','latex');            
       end 
       
       function [f,hP] = plotMeanVariable(obj,legendPlot,var)
            x = log2(obj.epsilonMeshSizeRatio(1,:));
            y = var';
            f = figure();
            hP{1} = plot(x,y,'+-');
            set(gca, 'XTickLabel',[])
            set(gca,'XTick',x,'XTickLabel',{'1' '2' '4' '8' '16' '32'})
            set(gca, 'Xdir', 'reverse');
            legObj = legend(legendPlot);            
            set(legObj,'Interpreter','latex','Location','best');
            xlabel('$\varepsilon/h$','interpreter','latex');                    
       end              
       
       function printPlot(obj,f,h,name)
           nameF = [name,obj.circleCase];
           outputName = strcat(obj.outputFolder,nameF);
           printer = plotPrinter(f,h);
           printer.print(outputName);             
       end              
        
    end
    
end