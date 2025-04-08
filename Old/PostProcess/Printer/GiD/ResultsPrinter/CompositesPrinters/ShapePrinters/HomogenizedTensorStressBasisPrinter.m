classdef HomogenizedTensorStressBasisPrinter < AbstractHomogenizedTensorPrinter

    methods (Access = public)
        
        function obj = HomogenizedTensorStressBasisPrinter(d)
            obj.simulationStr = 'HomogenizedTensorStressBasis';
            obj.computeNstre(d.ndim);
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function storeMicroProblemsFields(obj,d)
            microProblems = d.phyProblems{1};
            fields = microProblems.variables2printStressBasis;
            for istre = 1:obj.nstre
                di.fields = fields{istre};
                p = obj.printers{istre};
                p.storeFieldsToPrint(di);
                p.setStrVariablesNames('StressBasis');
                p.setStrVariablesMicroCase(istre);
            end
        end
    end
    
end