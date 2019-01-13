classdef ResultsPrinter < FilePrinter
    
    properties (Access = protected)
        ext = 'res'
        gaussDescriptor
        ngaus
        posgp
        ptype
        fields
        hasGaussData
    end
    
    properties (Access = protected, Abstract)
        simulationStr
    end
    
    methods (Access = public, Static)
        
        function rP = create(resultCase,d,dT)
            f = ResultsPrinterFactory();
            rP = f.create(resultCase,d,dT);
        end
        
    end
    
    methods (Access = public)
        
        function print(obj,iter,fields)
            obj.iter  = iter;
            obj.fields = fields;
            obj.createFileName(iter);
            obj.openFile();
            obj.printHeader();
            obj.printResults();
            obj.closeFile();
        end
        
        function printOnlyResults(obj,iter,fields)
            obj.openMode = 'a';
            obj.iter  = iter;
            obj.fields = fields;
            obj.createFileName(iter);
            obj.openFile();
            obj.printResults();
            obj.closeFile();
        end
        
        function f = getFieldName(obj)
            f = obj.fileName;
        end
        
        function setSimulationStr(obj,s)
            obj.simulationStr = s;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.storeDataBase(d)
            obj.createFileName(obj.iter);
            obj.openFile();
        end
        
        function storeDataBase(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                fieldName = fieldsNames{ifield};
                fieldValue = d.(fieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
        function d = createScalarGaussDataBase(obj,varargin)
            d = obj.createScalarDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
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
        
        function d = createVectorGaussDataBase(obj,varargin)
            d = obj.createVectorDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
        end
    end
    
    methods (Access = private)
        
        function d = addGaussDescriptor(obj,d)
            d.gaussDescriptor = obj.gaussDescriptor;
        end
        
    end
    
    methods (Abstract, Access = protected)
        printResults(obj)
        printHeader(obj)
    end
    
    
end
