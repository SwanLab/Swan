classdef ResultsPrinter < handle
    
    properties (Access = protected)
        fileID
        testName
        fileName
        nsteps
        etype
        ptype
        ngaus
        ndim
        posgp
        istep
        gaussDescriptor = 'Guass up?'
        
        ShapeNames %Refactor and delete
        optimizer
        printMode 
        hasGaussData

        
        fields
        openMode = 'w'
    end
    
    properties (Access = protected, Abstract)
        simulationStr
    end
    
    methods (Access = public)
        
        function print(obj,istep,fields)            
            obj.istep  = istep;
            obj.fields = fields;
            obj.createFileName();
            obj.openFile();
            obj.printHeader();            
            obj.printResults();
            obj.closeFile();
        end
        
        function printOnlyResults(obj,istep,fields)
            obj.openMode = 'a';
            obj.istep  = istep;
            obj.fields = fields;
            obj.createFileName();
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
        
        function setOpenMode(obj,s)
            obj.openMode = s;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.storeDataBase(d)
            obj.createFileName();
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
        
        function d = createScalarDataBase(obj,fieldValues,fieldName, fieldPosition)
            d.fileID = obj.fileID;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = obj.istep;
            d.fieldPosition = fieldPosition;
            d.simulationStr = obj.simulationStr;
        end
        
        function d = createScalarGaussDataBase(obj,varargin)
            d = obj.createScalarDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
        end
        
        function d = createVectorDataBase(obj,fieldValues,fieldName,fieldPosition,fieldComponentName)
            d = obj.createScalarDataBase(fieldValues,fieldName,fieldPosition);
            d.fieldComponentName = fieldComponentName;
        end
        
        function d = createVectorGaussDataBase(obj,varargin)
            d = obj.createVectorDataBase(varargin{:});
            d = obj.addGaussDescriptor(d);
        end
        
        function d = addGaussDescriptor(obj,d)
            d.gaussDescriptor = obj.gaussDescriptor;
        end
        
        function closeFile(obj)
            fclose(obj.fileID);
        end
        
        function createFileName(obj)
            iS = obj.istep;
            obj.fileName = fullfile('Output',obj.testName,strcat(obj.testName,num2str(iS),'.flavia.res'));
        end
        
        function openFile(obj)
            obj.fileID = fopen(obj.fileName,obj.openMode);
        end        
        
    end
    
    methods (Abstract, Access = protected)
        printResults(obj)
        printHeader(obj)
    end
    
    
end

