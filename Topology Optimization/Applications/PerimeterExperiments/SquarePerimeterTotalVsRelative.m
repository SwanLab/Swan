classdef SquarePerimeterTotalVsRelative < handle
    
    properties (Access = private)
        inputFile
        scale

        widthX
        widthY
        mesh
        levelSetCreator
        analyticTotalPerimeter
        analyticRelativePerimeter        
        regularizedRelativePerimeter
        regularizedTotalPerimeter
        geometricRelativePerimeter
        domainLengthX
        domainLengthY
        imesh
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
            for imesh = 1:numel(obj.inputFile)
                obj.imesh = imesh;
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
            obj.widthX = 0.50001;
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
            [connec,coord] = obj.loadSquareMeshParams();
            s.coord  = coord;
            s.connec = connec;
            obj.mesh = Mesh_Total(s);
        end
        
        function [connec,coord] = loadSquareMeshParams(obj)
            eval(obj.inputFile{obj.imesh});
            coord  = coord(:,2:3);
            connec = connec(:,2:end);
            obj.domainLengthX = max(coord(:,1)) - min(coord(:,1));
            obj.domainLengthY = max(coord(:,2)) - min(coord(:,2));            
        end
        
        function createLevelSet(obj)
            s.inputFile             = obj.inputFile{obj.imesh};
            s.mesh                  = obj.mesh;
            s.scale                 = obj.scale;
            s.plotting              = true;
            s.printing              = true;
            s.levelSetCreatorParams = obj.createLevelSetCreatorParams();
            lsCreator = LevelSetCreatorForPerimeter(s);
            obj.levelSetCreator = lsCreator;
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
            s.inputFile        = obj.inputFile{obj.imesh};
            s.mesh             = obj.mesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSetCreator.levelSet;
            s.outputFigureName = ['RelativeSmoothedRectangleMesh',num2str(obj.imesh)];
            s.plotting         = false;
            s.printing         = true;
            s.capturingImage   = true;
            s.isRobinTermAdded = false;
            s.perimeterType    = 'perimeterInterior';
            rPerimeter = RegularizedPerimeterComputer(s);
            rPerimeter.compute();
        end
        function tPerimeter = computeTotalPerimeters(obj)
            s.inputFile        = obj.inputFile{obj.imesh};
            s.mesh             = obj.mesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSetCreator.levelSet;
            s.outputFigureName = ['TotalSmoothedRectangleMesh',num2str(obj.imesh)];
            s.plotting         = false;
            s.printing         = true;
            s.capturingImage   = true;
            s.isRobinTermAdded = true;
            s.perimeterType    = 'perimeter';            
            tPerimeter = RegularizedPerimeterComputer(s);
            tPerimeter.compute();
        end        
        
        function gPerimeter = computeGeometricPerimeter(obj)
            s.designVariable = obj.levelSetCreator.levelSet;
            gPerimeter = GeometricPerimeterComputer(s);
            gPerimeter.compute();
        end
        
        function storePlotInfo(obj)
            h = obj.mesh.computeMeanCellSize;
            L = obj.domainLengthY;
            obj.xplot{obj.imesh} = obj.regularizedTotalPerimeter.epsilons/h;
            PerR = obj.regularizedRelativePerimeter.perimeters; 
            PerT = obj.regularizedTotalPerimeter.perimeters; 
            obj.relPerimeter{obj.imesh} = PerR;
           % obj.totPerimeter{obj.imesh} = (PerT - PerR)/2 + PerR;
             obj.totPerimeter{obj.imesh} = PerT;
            [n,d] = rat(h/L);
            obj.legendRelativePlot{obj.imesh} = ['$Per^R_\varepsilon \ (h/L = ',num2str(n),'/',num2str(d),')$'];
            obj.legendTotalPlot{obj.imesh} = ['$Per^T_\varepsilon \ (h/L = ',num2str(n),'/',num2str(d),')$'];
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
            printer.print(outputName);
        end
        
    end
    
end