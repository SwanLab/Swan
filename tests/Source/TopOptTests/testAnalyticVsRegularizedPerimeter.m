classdef testAnalyticVsRegularizedPerimeter < testShowingError ... & testTopOptComputation
    
    properties (Access = private)
        designVariable
        mesh
        levelSetCreator
        regularizedPerimeter
        analyticPerimeter
        radius
        inputFile
    end
    
    properties (Access = protected)
        testName = 'testAnalyticVsRegularizedPerimeter';
        tol
    end
    
    methods (Access = public)
        
        function obj = testAnalyticVsRegularizedPerimeter()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.computeRegularizedPerimeter();
        end
        
    end
    
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.designVariable = obj.topOpt.designVariable;
        end
        
        function computeError(obj)
            aP = obj.analyticPerimeter;
            nP = obj.regularizedPerimeter;
            obj.error = abs(aP - nP)/aP;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.tol = 5e-2;
            obj.radius = 0.31;
            obj.analyticPerimeter = 2*pi*obj.radius;   
            obj.inputFile = 'SquareMacroTriangle';
        end
        
         function createMesh(obj)
            [connec,coord] = obj.loadSquareMeshParams();
            s.coord  = coord;
            s.connec = connec;
            obj.mesh = Mesh_Total(s);
        end
        
        function [connec,coord] = loadSquareMeshParams(obj)
            eval(obj.inputFile);
            coord  = coord(:,2:3);
            connec = connec(:,2:end);
        end
        
        function createLevelSet(obj)
            s.inputFile    = obj.inputFile;
            s.mesh         = obj.mesh;
            s.scale        = 'MACRO';
            s.plotting              = false;
            s.printing              = false;
            s.levelSetCreatorParams = obj.createLevelSetCreatorParams();
            lsCreator = LevelSetCreatorForPerimeter(s);
            obj.designVariable = lsCreator.levelSet;
        end
        
        function s = createLevelSetCreatorParams(obj)
            ss.type = 'circleInclusion';
            domainLength = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            halfSide = domainLength/2;
            ss.fracRadius = (obj.radius/halfSide);            
            s = SettingsLevelSetCreator;
            s = s.create(ss); 
        end        
        
         function computeRegularizedPerimeter(obj)
            s = obj.createPerimeterParams();
            perimeterFunc = ShFunc_Perimeter(s);
            perimeterFunc.computeCostAndGradient();
            obj.regularizedPerimeter = perimeterFunc.value;             
        end     
        
          function s = createPerimeterParams(obj)
           sC.inputFile      = obj.inputFile;
           sC.mesh           = obj.mesh; 
           sC.designVariable = obj.designVariable;
           sC.epsilon        = obj.designVariable.mesh.computeMeanCellSize();
           sC.scale          = 'MACRO';
           fCreator = PerimeterParamsCreator(sC);        
           s = fCreator.perimeterParams;
        end     
                   
    end
    
end