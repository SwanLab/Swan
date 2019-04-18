classdef DesignVarMonitorFactory < handle
    
    properties (Access = private)
        monitor
        builder
        
        designVar
        optimizer
        dim
        showBC
    end
    
    
    methods (Access = public)
        
        function monitor = create(obj,cParams)
            shallDisplay = cParams.plotting;
            obj.init(cParams);
            if shallDisplay
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
            obj.optimizer = cParams.settings.optimizer;
            obj.dim = cParams.settings.pdim;
            obj.designVar = cParams.designVar;
            obj.showBC = cParams.settings.showBC;
        end
        
        function createMonitor(obj)
            switch obj.designVariable()
                case 'Density'
                    obj.monitor = DesignVarMonitor_Density(obj.designVar,obj.showBC);
                case 'LevelSet'
                    switch obj.dim
                        case '2D'
                            obj.monitor = DesignVarMonitor_LevelSet_2D(obj.designVar,obj.showBC);
                        case '3D'
                            obj.monitor = DesignVarMonitor_LevelSet_3D(obj.designVar,obj.showBC);
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
        
        function var = designVariable(obj)
            switch obj.optimizer
                case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                    var = 'LevelSet';
                case {'PROJECTED GRADIENT','MMA','IPOPT'}
                    var = 'Density';
            end
        end
        
        function returnNullMonitor(obj)
            obj.monitor = DesignVarMonitor_Null(obj.designVar,false);
        end
        
    end
    
end