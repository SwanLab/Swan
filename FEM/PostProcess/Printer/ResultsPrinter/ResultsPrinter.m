classdef ResultsPrinter < FilePrinter
    
    properties (Access = protected)
        ext = 'res'
        ptype
        fields
  end
    
    properties (Access = protected, Abstract)
        simulationStr
    end
    
    methods (Access = public, Static)
        
        function rP = create(resultCase,d)
            f = ResultsPrinterFactory();
            rP = f.create(resultCase,d);
        end
        
    end
    
    methods (Access = public)
        
        function print(obj,iter,d)
            obj.storeResultsInfo(d)
            obj.iter = iter;
            obj.createFileName(iter);
            obj.openFile();
            obj.printHeader();
            obj.printResults();
            obj.closeFile();
        end
        
        function printOnlyResults(obj,iter)
            obj.openMode = 'a';
            obj.iter  = iter;
            obj.createFileName(iter);
            obj.openFile();
            obj.printResults();
            obj.closeFile();
        end
        
        function setSimulationStr(obj,s)
            obj.simulationStr = s;
        end
        
    end
    
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.ptype = d.ptype;
            obj.etype = d.etype;
            obj.ndim  = d.ndim;
            obj.outFileName = d.outFileName;
            obj.resultsDir = d.resultsDir;            
        end
        
        function d = createScalarDataBase(obj,fieldValues,fieldName, fieldPosition)
            d.fileID = obj.fileID;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.iter = obj.iter;
            d.fieldPosition = fieldPosition;
            d.simulationStr = obj.simulationStr;
        end
        
        function d = createVectorDataBase(obj,fieldValues,fieldName,fieldPosition,fieldComponentName)
            d = obj.createScalarDataBase(fieldValues,fieldName,fieldPosition);
            d.fieldComponentName = fieldComponentName;
        end
        

    end
    
    methods (Abstract, Access = protected)
        printResults(obj)
        printHeader(obj)
    end
    
    methods (Access = public, Abstract)
       storeResultsInfo(obj)
    end
    
    
end
