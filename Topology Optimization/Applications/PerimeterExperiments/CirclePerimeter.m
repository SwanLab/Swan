classdef CirclePerimeter < handle
    
    properties (Access = private)
        inputFile
        scale
        radius
        
        mesh
        levelSet
        analyticPerimeter
        regularizedPerimeter
        geometricPerimeter
        domainLength
        imesh
        xplot
        yplot
        legendPlot
        nameCase
        outputFolder
    end
    
    methods (Access = public)
        
        function obj = CirclePerimeter()
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
            obj.radius = 0.21;
            obj.nameCase = 'CirclePerimeterExperiment';
            obj.inputFile = {'SquareMacroTriangle'; ...
                             'SquareMacroTriangleFine';...
                             'SquareMacroTriangleFineFine'};
            obj.scale     = 'MACRO';
            obj.analyticPerimeter    = 2*pi*(obj.radius);
            obj.outputFolder = '/home/alex/Dropbox/Perimeter/';
        end
        
        function createMesh(obj)
            s.inputFile = obj.inputFile{obj.imesh};
            meshCreator = MeshCreatorFromInputFile(s);
            meshCreator.createMesh();
            obj.mesh         = meshCreator.mesh;
            obj.domainLength = meshCreator.domainLength;
        end
        
        function createLevelSet(obj)
            s.inputFile             = obj.inputFile{obj.imesh};
            s.mesh                  = obj.mesh;
            s.scale                 = obj.scale;
            s.plotting              = true;
            s.printing              = true;
            s.levelSetCreatorParams = obj.createLevelSetCreatorParams();
            lsCreator = LevelSetCreatorForPerimeter(s);
            obj.levelSet = lsCreator.levelSet;
        end
        
        function s = createLevelSetCreatorParams(obj)
            ss.type = 'circleInclusion';
            halfSide = obj.domainLength/2;
            ss.fracRadius = (obj.radius/halfSide);            
            s = SettingsLevelSetCreator;
            s = s.create(ss); 
        end
        
        function computePerimeters(obj)
            obj.regularizedPerimeter = obj.computeRegularizedPerimeters();
            obj.geometricPerimeter   = obj.computeGeometricPerimeter();
        end
        
        function rPerimeter = computeRegularizedPerimeters(obj)
            s.inputFile        = obj.inputFile{obj.imesh};
            s.mesh             = obj.mesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSet;
            s.outputFigureName = ['SmoothedCircleMesh',num2str(obj.imesh)];
            s.plotting         = true;
            s.printing         = true;
            s.capturingImage   = true;
            s.isRobinTermAdded = false;
            s.perimeterType    = 'perimeterInterior';            
            rPerimeter = RegularizedPerimeterComputer(s);
            rPerimeter.compute();
        end
        
        function gPerimeter = computeGeometricPerimeter(obj)
            s.designVariable = obj.levelSet;
            gPerimeter = GeometricPerimeterComputer(s);
            gPerimeter.compute();
        end
        
        function storePlotInfo(obj)
            h = obj.mesh.computeMeanCellSize;
            L = obj.domainLength;
            obj.xplot{obj.imesh} = obj.regularizedPerimeter.epsilons/h;
            obj.yplot{obj.imesh} = obj.regularizedPerimeter.perimeters;
            [n,d] = rat(h/L);
            obj.legendPlot{obj.imesh} = ['$Per^R_\varepsilon \ (h/L = ',num2str(n),'/',num2str(d),')$'];
        end
        
        function plotPerimeterValues(obj)
            f = figure();
            hold on
            nmesh = numel(obj.inputFile);
            p = cell(nmesh+1,1);
            leg = cell(nmesh+1,1);
            leg{1} = '$2\pi r$';
            p{1} = plot(log2(obj.xplot{1}),obj.analyticPerimeter*ones(size(obj.xplot{1})));            
            for im = 1:nmesh
                x = log2(obj.xplot{im});
                y = obj.yplot{im};
                hold on
                p{im+1}   = plot(x,y,'-+');
                leg{im+1} = obj.legendPlot{im};
                set(gca, 'XTickLabel',[])                      
                set(gca,'XTick',x,'XTickLabel',{'1' '2' '4' '8' '16' '32'})                
            end
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','southeast');
            xlabel('$\varepsilon/h$','interpreter','latex');
            ylim([0 1.1*obj.analyticPerimeter])
            set(gca, 'Xdir', 'reverse');
            obj.printPlot(f,p);
        end
        
        function printPlot(obj,f,h)
            outputName = strcat(obj.outputFolder,obj.nameCase);
            printer = plotPrinter(f,h);
            printer.print(outputName);
        end
        
    end
    
end