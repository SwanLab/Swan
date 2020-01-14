classdef DensityPlotterForPerimeter < handle
    
    properties (Access = private)
        inputFile
        mesh
        scale
        density
        densityVariable
        plotter
    end
    
    methods (Access = public)
        
        function obj = DensityPlotterForPerimeter(cParams)
            obj.init(cParams)
        end
        
        function plot(obj)  
            obj.createDesignVariable();
            obj.createPlotter();
            obj.plotter.refresh();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.inputFile  = cParams.inputFile;
            obj.scale      = cParams.scale;
            obj.density    = cParams.density;
        end
        
        function createDesignVariable(obj)
            s = SettingsDesignVariable();
            s.mesh                    = obj.mesh;
            s.type                    = 'Density';
            s.initialCase             = 'full';
            s.levelSetCreatorSettings = obj.createLevelSetParams();
            s.scalarProductSettings   = obj.createScalarProductParams();
            s.femData                 = obj.createFemContainerData();
            obj.densityVariable       = DesignVariable.create(s);            
            obj.densityVariable.value = obj.density;
            obj.densityVariable.rho   = obj.density;                        
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
        
        function s = createScalarProductParams(obj)
            s.scalarProductSettings.femSettings = [];
            s.epsilon = [];     
        end
        
        function s = createLevelSetParams(obj)
            ss.type = 'full';
            s = SettingsLevelSetCreator;
            s = s.create(ss);
        end        
        
        function createPlotter(obj)
            sD.shallDisplay   = true;
            sD.showBC         = false;
            sD.bc             = false;
            sD.designVariable = obj.densityVariable;
            sD.optimizerName  = '';
            sD.dim            = '2D';
            sD.scale          = obj.scale;
            f = DesignVarMonitorFactory;
            obj.plotter = f.create(sD);                        
        end
        
    end    
    
end