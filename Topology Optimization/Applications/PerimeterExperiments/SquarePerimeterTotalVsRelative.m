classdef SquarePerimeterTotalVsRelative < handle
    
    properties (Access = private)
        inputFile
        scale

        widthX
        widthY
        backgroundMesh
        boundaryMesh
        levelSet
        analyticTotalPerimeter
        analyticRelativePerimeter        
        regularizedRelativePerimeter
        regularizedTotalPerimeter
        geometricRelativePerimeter
        domainLengthX
        domainLengthY
        iMesh
        xplot
        totPerimeter
        relPerimeter
        legendRelativePlot
        legendTotalPlot
        nameCase
        outputFolder
    end
    
    methods (Access = public)
        
        function obj = SquarePerimeterTotalVsRelative()
            obj.init();
            for im = 1:numel(obj.inputFile)
                obj.iMesh = im;
                obj.createMesh();
                obj.createLevelSet();
                obj.computePerimeters();
                obj.storePlotInfo();
            end
            obj.plotPerimeterValues();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.widthX = 0.52;
            obj.widthY = 1;
            %obj.inputFile = 'DoubleSquareMacroTriangle';
            obj.nameCase = 'SquareTotalVsRelativePerimeterExperiment';
            obj.inputFile = {'RectangleMacroTriangle';...
                             'RectangleMacroTriangleFine';...
                             'RectangleMacroTriangleFineFine'};
            obj.scale     = 'MACRO';
            %obj.analyticPerimeter    = 0.7*4;
            obj.outputFolder = '/home/alex/Dropbox/Perimeter/';
        end
        
        function createMesh(obj)
            s.inputFile = obj.inputFile{obj.iMesh};
            s.isBackgroundMeshRectangularBox = true;
            mCreator = BackgroundAndBoundaryMeshCreatorFromInputFile(s);
            obj.backgroundMesh = mCreator.backgroundMesh;
            obj.boundaryMesh = mCreator.boundaryMesh;
        end
        
        function loadSquareMeshParams(obj)
            coord  = obj.backgroundMesh.coord;
            obj.domainLengthX = max(coord(:,1)) - min(coord(:,1));
            obj.domainLengthY = max(coord(:,2)) - min(coord(:,2));            
        end
        
        function createLevelSet(obj)
            s = obj.createLevelSetCreatorParams();
            s.coord      = obj.backgroundMesh.coord;
            s.ndim       = obj.backgroundMesh.ndim;
            lsCreator = LevelSetCreator.create(s);
            obj.levelSet = lsCreator.getValue();            
        end
        
        function s = createLevelSetCreatorParams(obj)
            ss.widthH = obj.widthX;
            ss.widthV = 1.1*obj.widthY;                        
            ss.type = 'rectangle';
            s = SettingsLevelSetCreator;
            s = s.create(ss); 
        end
        
        function computePerimeters(obj)
            obj.analyticTotalPerimeter    = 2*obj.widthX*obj.domainLengthX + 2*obj.widthY*obj.domainLengthY;
            obj.analyticRelativePerimeter = 2*obj.widthY*obj.domainLengthY;            
            obj.regularizedRelativePerimeter = obj.computeRelativePerimeters();
            obj.regularizedTotalPerimeter    = obj.computeTotalPerimeters();            
            obj.geometricRelativePerimeter   = obj.computeGeometricPerimeter();
        end
        
        function rPerimeter = computeRelativePerimeters(obj)
            s.inputFile        = obj.inputFile{obj.iMesh};
            s.backgroundMesh   = obj.backgroundMesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSet;
            s.outputFigureName = ['RelativeSmoothedRectangleMeshColor',num2str(obj.iMesh)];
            s.plotting         = false;
            s.printing         = true;
            s.capturingImage   = true;
            s.isRobinTermAdded = false;
            s.perimeterType    = 'perimeterInterior';
            rPerimeter = RegularizedPerimeterComputer(s);
            rPerimeter.compute();
        end
        
        function tPerimeter = computeTotalPerimeters(obj)
            s.inputFile        = obj.inputFile{obj.iMesh};
            s.backgroundMesh   = obj.backgroundMesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSet;
            s.outputFigureName = ['TotalSmoothedRectangleMeshColor',num2str(obj.iMesh)];
            s.plotting         = false;
            s.printing         = true;
            s.capturingImage   = true;
            s.isRobinTermAdded = true;
            s.perimeterType    = 'perimeter';            
            tPerimeter = RegularizedPerimeterComputer(s);
            tPerimeter.compute();
        end        
        
        function gPerimeter = computeGeometricPerimeter(obj)
            s.unfittedMesh = obj.createUnfittedMesh();            
            gPerimeter = GeometricPerimeterComputer(s);
            gPerimeter.compute();
        end
        
        function uMesh = createUnfittedMesh(obj)
            sM.backgroundMesh = obj.backgroundMesh;
            sM.boundaryMesh   = obj.boundaryMesh;            
            uMesh = UnfittedMesh(sM);
            uMesh.compute(obj.levelSet);
        end
        
        function storePlotInfo(obj)
            obj.loadSquareMeshParams();            
            h = obj.backgroundMesh.computeMeanCellSize;
            L = obj.domainLengthY;
            obj.xplot{obj.iMesh} = obj.regularizedTotalPerimeter.epsilons/h;
            PerR = obj.regularizedRelativePerimeter.perimeters; 
            PerT = obj.regularizedTotalPerimeter.perimeters; 
            obj.relPerimeter{obj.iMesh} = PerR;
           % obj.totPerimeter{obj.iMesh} = (PerT - PerR)/2 + PerR;
             obj.totPerimeter{obj.iMesh} = PerT;
            [n,d] = rat(h/L);
            obj.legendRelativePlot{obj.iMesh} = ['$Per^R_\varepsilon \ (h/L = ',num2str(n),'/',num2str(d),')$'];
            obj.legendTotalPlot{obj.iMesh} = ['$Per^T_\varepsilon \ (h/L = ',num2str(n),'/',num2str(d),')$'];
        end
        
        function plotPerimeterValues(obj)
            obj.plotPerimeterTotal();
            obj.plotPerimeterRelative();            
        end        
        
        function plotPerimeterRelative(obj)
            s.analyticPerimeter       = obj.analyticRelativePerimeter;
            s.analyticPerimeterLegend = '$2L$';
            s.titlePlot               = 'RelativePer';
            s.perimeterValues         = obj.relPerimeter;
            s.legendPlot              = obj.legendRelativePlot;
            obj.plotPerimeter(s)                 
        end
        
        function plotPerimeterTotal(obj)
            s.analyticPerimeter       = obj.analyticTotalPerimeter;
            s.analyticPerimeterLegend = '$4L$';
            s.titlePlot               = 'TotalPer';
            s.perimeterValues         = obj.totPerimeter;
            s.legendPlot              = obj.legendTotalPlot;
            obj.plotPerimeter(s)                                   
        end                
        
        function plotPerimeter(obj,cParams)
            aPerimeter        = cParams.analyticPerimeter*ones(size(obj.xplot{1}));
            legendAnalyticPer = cParams.analyticPerimeterLegend;
            titlePlot         = cParams.titlePlot;
            perimeterValues   = cParams.perimeterValues;
            legends           = cParams.legendPlot;
            f = figure();
            hold on
            nmesh = numel(obj.inputFile);
            p = cell(nmesh+1,1);
            leg = cell(nmesh+1,1);
            leg{1} = legendAnalyticPer;
            p{1} = plot(log2(obj.xplot{1}),aPerimeter);     
            for im = 1:nmesh
                x = log2(obj.xplot{im});
                y = perimeterValues{im};
                hold on
                p{im+1}   = plot(x,y,'-+');
                leg{im+1} = legends{im};
                set(gca, 'XTickLabel',[])                      
                set(gca,'XTick',x,'XTickLabel',{'1' '2' '4' '8' '16' '32'})                
            end
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','southeast');
            xlabel('$\varepsilon/h$','interpreter','latex');
           % ylim([0 1.1*max(aPerimeter(:))])
            set(gca, 'Xdir', 'reverse');
            obj.printPlot(f,p,titlePlot);                  
        end

        
        function printPlot(obj,f,h,perCase)
            outputName = strcat(obj.outputFolder,perCase,obj.nameCase);
            printer = plotPrinter(f,h);
           % printer.print(outputName);
        end
        
    end
    
end