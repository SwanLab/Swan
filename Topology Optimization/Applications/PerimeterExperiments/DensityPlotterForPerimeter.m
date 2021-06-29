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
            s.mesh      = obj.mesh;
            s.inputFile = obj.inputFile;
            s.scale     = obj.scale;
            s.type      = 'Density';
            d = DesignVariableCreatorSettings(s);
            s = d.create();
            obj.densityVariable       = DesignVariable.create(s);            
            obj.densityVariable.update(obj.density);
            obj.densityVariable.rho   = obj.density;                        
        end
        
    
        
        function createPlotter(obj)
            sD.shallDisplay   = true;
            sD.showBC         = false;
            sD.bc             = false;
            sD.designVariable = obj.densityVariable;
            sD.optimizerNames = '';
            sD.dim            = '2D';
            sD.scale          = obj.scale;
            sD.mesh           = Mesh_Total(obj.mesh);
            f = DesignVarMonitorFactory;
            obj.plotter = f.create(sD);                        
        end
        
    end    
    
end