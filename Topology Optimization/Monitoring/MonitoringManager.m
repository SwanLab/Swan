classdef MonitoringManager < handle
    
    properties (Access = public)
        monitoring
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MonitoringManager(cParams)
            switch cParams.outputFunction.type
                case 'Academic'
                    obj.monitoring = AcademicMonitoring(cParams);
                case 'Topology'
                    obj.monitoring = TopologyMonitoring(cParams);
                case 'Truss Structure'
                    obj.monitoring = TrussStructureMonitoring(cParams);
            end
        end
        
    end
    
end