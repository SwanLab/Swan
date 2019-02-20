classdef DesignVarMonitorFactory < handle
    
    properties (Access = private)
        monitor
        builder
        
        mesh
        optimizer
        dim
        showBC
    end
    
    
    methods (Access = public)
        
        function monitor = create(obj,shallDisplay,settings,mesh)
            obj.init(settings,mesh);
            if shallDisplay
                obj.createBuilder();
                obj.createMonitor();
                obj.build();
            else
                obj.returnNullMonitor(mesh);
            end
            
            monitor = obj.monitor;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,settings,mesh)
            obj.optimizer = settings.optimizer;
            obj.dim = settings.pdim;
            obj.mesh = mesh;
            obj.showBC = settings.showBC;
        end
        
        function createMonitor(obj)
            switch obj.designVariable()
                case 'Density'
                    obj.monitor = DesignVarMonitor_Density(obj.mesh,obj.showBC);
                case 'LevelSet'
                    switch obj.dim
                        case '2D'
                            obj.monitor = DesignVarMonitor_LevelSet_2D(obj.mesh,obj.showBC);
                        case '3D'
                            obj.monitor = DesignVarMonitor_LevelSet_3D(obj.mesh,obj.showBC);
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
        
        function returnNullMonitor(obj,mesh)
            obj.monitor = DesignVarMonitor_Null(mesh);
        end
        
    end
    
end