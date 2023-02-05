classdef DesignVarMonitorFactory < handle
    
    properties (Access = private)
        monitor
        builder
        
        scale
        designVariable
        sectionVariables
        optimizer
        dim
        showBC
        bc
        shallDisplay
        mesh
    end
    
    
    methods (Access = public)
        
        function monitor = create(obj,cParams)
            obj.init(cParams);
            if obj.shallDisplay
                obj.createBuilder();
                obj.createMonitor();
                obj.build();
            else
                obj.returnNullMonitor();
            end
            
            monitor = obj.monitor;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.shallDisplay   = cParams.shallDisplay;
            obj.optimizer      = cParams.optimizerNames;
            obj.dim            = cParams.dim;
            obj.designVariable = cParams.designVariable;
            obj.sectionVariables = cParams.sectionVariables;
            obj.showBC         = cParams.showBC;
            obj.bc             = cParams.bc;
            obj.scale          = cParams.scale;
            obj.mesh           = cParams.mesh;
        end
        
        function createMonitor(obj)
            mS = obj.getMonitorSettings();
            switch obj.designVariable.type
                case {'Density','MicroParams','DensityEigModes'}
                    obj.monitor = DesignVarMonitor_Density(mS);
                case {'LevelSet','LevelSetEigModes'}
                    switch obj.dim
                        case '2D'
                            obj.monitor = DesignVarMonitor_LevelSet_2D(mS);
                        case '3D'
                            obj.monitor = DesignVarMonitor_LevelSet_3D(mS);
                    end
                case 'AreaColumn'
                    obj.monitor = DesignVarMonitor_AreaColumn(mS);
                case 'RadiusColumn'
                    obj.monitor = DesignVarMonitor_RadiusColumn(mS);
                case 'SquareColumn'
                    obj.monitor = DesignVarMonitor_SquareColumn(mS);
                case 'RectangularColumn'
                    obj.monitor = DesignVarMonitor_RectangularColumn(mS);
                case 'HoleColumn'
                    obj.monitor = DesignVarMonitor_HoleColumn(mS);
                otherwise
                    error('Invalid Design Variable')
            end
        end
        
        function createBuilder(obj)
            obj.builder = Factory_BuilderDesignVarMonitor().create(obj.mesh.ndim);
        end
        
        function build(obj)
            obj.builder.build(obj.monitor);
        end
        
        function returnNullMonitor(obj)
            mS = obj.getMonitorSettings();
            obj.monitor  = DesignVarMonitor_Null(mS);
        end
        
        function s = getMonitorSettings(obj)
            s.designVar = obj.designVariable;
            s.sectionVariables = obj.sectionVariables;
            s.showBC    = obj.showBC;
            s.bc        = obj.bc;
            s.dim       = obj.dim;
            s.scale     = obj.scale;
            s.mesh      = obj.mesh;
        end
        
    end
    
end