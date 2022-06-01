classdef VectorFieldPrinter < ResultsPrinter ...
    
    properties (Access = protected)
        simulationStr = 'VectorFieldResults';
        hasGaussData = true;
    end
    
    properties (Access = private)
        dispStr   = 'Displacements';
        displCompStr  = 'U';
        dV
        nameCase
    end
    
    methods (Access = public)
        
        function obj = VectorFieldPrinter(d)
            obj.init(d);
            obj.nameCase = d.name;
            obj.setStrVariablesNames()
        end
        
        function printResults(obj,iter,fileID)
            obj.createDataBases(obj.fields,iter,fileID);
            VectorNodalPrinter(obj.dV);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.fields;
        end
        
    end
    
    methods (Access = private)
        
        function setStrVariablesNames(obj)
            n = obj.nameCase;
            obj.dispStr   = [obj.dispStr,n];
        end
        
        function createDataBases(obj,fields,iter,fileID)
            f = fields;
            obj.dV   = obj.createVectorDataBase(iter,fileID,f.u,obj.dispStr,'OnNodes',obj.displCompStr);
        end
        
    end
    
end

