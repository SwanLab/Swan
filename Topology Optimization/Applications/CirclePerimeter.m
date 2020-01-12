classdef CirclePerimeter < handle
    
    properties (Access = private)
        inputFile
        cParamsPerimeter
        cParamsFem
        cParamsFilter
        mesh
        designVariable
        homogenizedVarComputer
        targetParameters
        regularizedPerimeter
        radius
        unfittedMesh
        epsilon
    end
    
    methods (Access = public)
        
        function obj = CirclePerimeter()
            obj.init()
            obj.createPerimeterParams();
            obj.createPerimeterFunction();
            obj.createUnfittedMesh();
            obj.plotLevelSet();
            obj.printLevelSet();
            obj.computePerimeters();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.radius = 0.47;
            obj.inputFile = 'SquareMacroTriangleFineFine';
        end
        
        function createMesh(obj)
            [connec,coord] = loadSquareMeshParams(obj);
            s.coord  = coord;
            s.connec = connec;
            obj.mesh = Mesh_Total(s);
        end
        
        function [connec,coord] = loadSquareMeshParams(obj)
            eval(obj.inputFile);
            coord  = coord(:,2:3);
            connec = connec(:,2:end);
        end
        
        function createDesignVariable(obj)
            s = SettingsDesignVariable();
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.initialCase = 'full';
            s.levelSetCreatorSettings = obj.createLevelSetParams();
            s.scalarProductSettings   = obj.createScalarProductParams();
            s.femData = obj.createFemContainerData();
            s.scalarProductSettings.epsilon = 0.01;
            obj.designVariable = DesignVariable.create(s);
        end
        
        function s = createLevelSetParams(obj)
            ss.type = 'circleInclusion';
            ss.fracRadius = (obj.radius/0.5);
            s = SettingsLevelSetCreator;
            s = s.create(ss);
        end
        
        function s = createScalarProductParams(obj)
            s.scalarProductSettings.femSettings = obj.inputFile;
            s.scale = 'MACRO';
        end
        
        function s = createFemContainerData(obj)
            s = FemDataContainer;
            s.fileName = obj.inputFile;
            s.scale    = 'MACRO';
            s.pdim     = '2D';
            s.ptype    = 'ELASTIC';
            s.nelem    = size(obj.mesh.connec,1);
            s.bc       = [];
            s.coord    = obj.mesh.coord;
            s.connec   = obj.mesh.connec;
        end
        
        function createFemParams(obj)
            s.fileName = obj.inputFile;
            s.scale    = 'MACRO';
            s.mesh     = obj.mesh;
            s.isRobinTermAdded = false;
            obj.cParamsFem = s;
        end
        
        function createFilterParams(obj)
            s = SettingsFilter();
            s.filterType =  'PDE';
            s.domainType =  'INTERIOR';
            s.designVar  =  obj.designVariable;
            s.quadratureOrder =  'LINEAR';
            s.femSettings = obj.cParamsFem;
            obj.cParamsFilter = s;
        end
        
        function createHomogenizedVarComputer(obj)
            s = SettingsHomogenizedVarComputer;
            cParams = s;
            obj.homogenizedVarComputer = HomogenizedVarComputer(cParams);
        end
        
        function s = createTargetParamters(obj)
            s = TargetParameters();
            s.epsilon_perimeter = obj.epsilon;
            s.epsilon = obj.epsilon;
        end
        
        function createPerimeterParams(obj)
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFemParams();
            obj.createFilterParams();                        
            s = SettingsShapeFunctional();
            s.filterParams = obj.cParamsFilter;
            s.femSettings  = obj.cParamsFem;
            s.homogVarComputer = obj.homogenizedVarComputer;
            s.designVariable   = obj.designVariable;
            s.targetParameters = obj.createTargetParamters();
            s.type = 'perimeterInterior';
            obj.cParamsPerimeter = s;
        end
        
        function createPerimeterFunction(obj)
            s = obj.cParamsPerimeter;
            shFunc = ShFunc_Perimeter(s);
            obj.regularizedPerimeter = shFunc;
        end
        
        function createUnfittedMesh(obj)
            cParams = SettingsMeshUnfitted('INTERIOR',obj.mesh);
            uMesh = UnfittedMesh(cParams);
            uMesh.compute(obj.designVariable.value);
            obj.unfittedMesh = uMesh;
        end
        
        function plotLevelSet(obj)
            obj.unfittedMesh.plot();
        end
        
        function printLevelSet(obj)
            sP.mesh = obj.designVariable.mesh;
            sP.outName = 'PerimeterExperimentLevelSet';
            sP.pdim  = '2D';
            sP.ptype = 'TRIANGLE';
            p = PostProcessDataBaseCreator(sP);
            s = p.getValue();
            type = 'LevelSet';
            printer = Postprocess(type,s);
            d.x = obj.designVariable.value;
            iter = 0;
            printer.print(iter,d);                
        end
        
        function Per = computeRegularizedPerimeters(obj)
            nepsilon = 5;
            epsmin = obj.designVariable.mesh.computeMeanCellSize;
            epsmax = 10*epsmin;
            eps = linspace(epsmax,epsmin,nepsilon);
            Per = zeros(nepsilon,1);
            for iepsilon = 1:nepsilon
               obj.epsilon = eps(iepsilon);
               Per(iepsilon) = obj.computeRegularizedPerimeter();
               obj.plotDensity();
               obj.printDensity(iepsilon);
            end   
        end
        
        function plotDensity(obj)
            s = SettingsDesignVariable();
            s.mesh = obj.mesh;
            s.type = 'Density';
            s.initialCase = 'full';
            s.levelSetCreatorSettings = obj.createLevelSetParams();
            s.scalarProductSettings   = obj.createScalarProductParams();
            s.femData = obj.createFemContainerData();
            s.scalarProductSettings.epsilon = 0.01;
            dVariable = DesignVariable.create(s);
            
            dVariable.value = obj.regularizedPerimeter.regularizedDensity;
            dVariable.rho = obj.regularizedPerimeter.regularizedDensity;
            
            sD.shallDisplay   = true;
            sD.showBC         = false;
            sD.bc             = false;
            sD.designVariable = dVariable;
            sD.optimizerName  = '';
            sD.dim            = '2D'; 
            sD.scale          = 'MACRO';
            f = DesignVarMonitorFactory;
            designVarMonitor = f.create(sD);            
            designVarMonitor.refresh();            
        end
        
        function printDensity(obj,iepsilon)
            sP.mesh = obj.designVariable.mesh;
            sP.outName = 'PerimeterExperiment';
            sP.pdim  = '2D';
            sP.ptype = 'TRIANGLE';
            p = PostProcessDataBaseCreator(sP);
            s = p.getValue();
            type = 'Density';
            printer = Postprocess(type,s);
            d.x = obj.regularizedPerimeter.regularizedDensity;
            iter = iepsilon;
            printer.print(iter,d);            
        end
        
        function Per = computeRegularizedPerimeter(obj)
            obj.createPerimeterParams();
            obj.createPerimeterFunction();      
            obj.regularizedPerimeter.computeCostAndGradient();
            Per = obj.regularizedPerimeter.value;            
        end
        
        function computePerimeters(obj)
            PerAnalytic = 2*pi*obj.radius
            PerNumericReguarized = obj.computeRegularizedPerimeters()'
            PerNumericGeometric  = obj.computeGeometricPerimeter()
        end
        
        function Per = computeGeometricPerimeter(obj)
            meshBackground = obj.designVariable.mesh;
            interpolation = Interpolation.create(meshBackground,'LINEAR');
            s.unfittedType = 'BOUNDARY';
            s.meshBackground = meshBackground;
            s.interpolationBackground = interpolation;
            s.includeBoxContour = false;
            cParams = SettingsMeshUnfitted(s);
            bMesh = UnfittedMesh(cParams);
            bMesh.compute(obj.designVariable.value);
            
            s.mesh = bMesh;
            s.type = 'COMPOSITE';
            sc.mesh = bMesh.innerCutMesh;
            sc.type = 'CutMesh';
            s.compositeParams{1} = sc;
            integrator = Integrator.create(s);
            npnod = bMesh.meshBackground.npnod;
            f = ones(npnod,1);
            int = integrator.integrateAndSum(f);
            Per = sum(int);
        end
        
        function plotUnfittedMesh(obj)
            obj.unfittedMesh.plot();
        end
        
    end
    
end