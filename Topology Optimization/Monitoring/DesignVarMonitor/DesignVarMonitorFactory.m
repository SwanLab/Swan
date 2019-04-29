classdef DesignVarMonitorFactory < handle
    
    properties (Access = private)
        monitor
        builder
        
        designVar
        optimizer
        dim
        showBC
        shallDisplay
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
            obj.shallDisplay = cParams.shallDisplay;
            obj.optimizer    = cParams.optimizerName;
            obj.dim          = cParams.dim;
            obj.designVar    = cParams.designVar;
            obj.showBC       = cParams.showBC;
        end
        
        function createMonitor(obj)
            mS.designVar = obj.designVar;
            mS.showBC    = obj.showBC;
            switch obj.designVar.type
                case 'Density'
                    obj.monitor = DesignVarMonitor_Density(mS);
                case 'LevelSet'
                    switch obj.dim
                        case '2D'
                            obj.monitor = DesignVarMonitor_LevelSet_2D(mS);
                        case '3D'
                            obj.monitor = DesignVarMonitor_LevelSet_3D(mS);
                    end
                otherwise
                    error('Invalid Design Variable')
            end
        end
        
        function createBuilder(obj)
            obj.builder = Factory_BuilderDesignVarMonitor().create(obj.dim);
        end
        
        function build(obj)
            obj.builder.build(obj.monitor);
        end
        
        function returnNullMonitor(obj)
            mS.designVar = obj.designVar;
            mS.showBC    = false;
            obj.monitor  = DesignVarMonitor_Null(mS);
        end
        
    end
    
end