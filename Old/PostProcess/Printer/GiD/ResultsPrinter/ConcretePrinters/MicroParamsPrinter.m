classdef MicroParamsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'NodalMicroParams';
        hasGaussData = false;
    end

    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MicroParamsPrinter(d)
            obj.init(d);
        end
        
        function printResults(obj,iter,fileID)
            obj.createDensityDataBase(iter,fileID);
            obj.createM1DataBase(iter,fileID);
            obj.createM2DataBase(iter,fileID);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            obj.fields = d.fields;
        end
        
    end
    
    methods (Access = private)
       
        function createM1DataBase(obj,iter,fileID)
            f = obj.fields{1};
            dS = obj.createScalarDataBase(iter,fileID,f,'M1','OnNodes');
            ScalarNodalPrinter(dS);
        end
        
        function createM2DataBase(obj,iter,fileID)
            f = obj.fields{2};
            dS = obj.createScalarDataBase(iter,fileID,f,'M2','OnNodes');
            ScalarNodalPrinter(dS);
        end
        
        function createDensityDataBase(obj,iter,fileID)
            f = obj.fields{3};
            dS = obj.createScalarDataBase(iter,fileID,f,'Density','OnNodes');
            ScalarNodalPrinter(dS);
        end
        
    end
    
end