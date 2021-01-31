classdef HomogenizedTensorPrinter < AbstractHomogenizedTensorPrinter
    
    methods (Access = public)
        
        function obj = HomogenizedTensorPrinter(d)
            obj.simulationStr = 'HomogenizedTensor';            
            obj.computeNstre(d.ndim);
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
               
        function storeMicroProblemsFields(obj,d)
            microProblems = d.phyProblems{1};
            fields = microProblems.variables2print;
            for istre = 1:obj.nstre
                di.fields = fields{istre};
                p = obj.printers{istre};
                p.storeFieldsToPrint(di);
                p.setStrVariablesMicroCase(istre)                                
            end                       
        end               
        
    end    
    
end