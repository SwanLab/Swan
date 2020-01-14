classdef LevelSetCreatorForPerimeter < handle
    
    properties (Access = public)
        levelSet
    end
    
    properties (Access = private)        
        inputFile
        mesh
        scale
        domainLength
        plotting
        printing
        levelSetCreatorParams
    end
    
    methods (Access = public)
        
        function obj = LevelSetCreatorForPerimeter(cParams)
            obj.init(cParams);
            obj.createLevelSet();
            obj.plotLevelSet();
            obj.printLevelSet();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.inputFile             = cParams.inputFile;
            obj.mesh                  = cParams.mesh;
            obj.scale                 = cParams.scale;
            obj.plotting              = cParams.plotting;
            obj.printing              = cParams.printing;
            obj.levelSetCreatorParams = cParams.levelSetCreatorParams;
        end
        
        function createLevelSet(obj)
            s = obj.createLevelSetParams();
            obj.levelSet = DesignVariable.create(s);            
        end
        
        function s = createLevelSetParams(obj)
            s = SettingsDesignVariable();
            s.mesh                    = obj.mesh;
            s.type                    = 'LevelSet';
            s.initialCase             = 'full';
            s.levelSetCreatorSettings = obj.levelSetCreatorParams;
            s.scalarProductSettings   = obj.createScalarProductParams();
            s.femData                 = obj.createFemContainerData();            
        end
        
        function s = createScalarProductParams(obj)            
            s.scalarProductSettings.femSettings = obj.inputFile;
            s.scale = obj.scale;
            s.epsilon = 0.01;            
        end
        
        function s = createFemContainerData(obj)
            s = FemDataContainer;
            s.fileName = obj.inputFile;
            s.scale    = obj.scale;
            s.pdim     = '2D';
            s.ptype    = 'ELASTIC';
            s.nelem    = size(obj.mesh.connec,1);
            s.bc       = [];
            s.coord    = obj.mesh.coord;
            s.connec   = obj.mesh.connec;
        end
        
        function plotLevelSet(obj)
            uMesh = obj.createUnfittedMesh();
            if obj.plotting
                uMesh.plot();
            end
        end                
        
       function mesh = createUnfittedMesh(obj)
            cParams = SettingsMeshUnfitted('INTERIOR',obj.mesh);
            mesh = UnfittedMesh(cParams);
            mesh.compute(obj.levelSet.value);           
        end
        
        function printLevelSet(obj)
            if obj.printing
            s = obj.createLevelSetPrinterParams();
            type = 'LevelSet';
            printer = Postprocess(type,s);
            d.x = obj.levelSet.value;
            iter = 0;
            printer.print(iter,d);                
            end
        end       
        
        function s = createLevelSetPrinterParams(obj)
            sP.mesh = obj.levelSet.mesh;
            sP.outName = 'PerimeterExperimentLevelSet';
            sP.pdim  = '2D';
            sP.ptype = 'TRIANGLE';
            p = PostProcessDataBaseCreator(sP);    
            s = p.getValue();            
        end
        
    end
    
end