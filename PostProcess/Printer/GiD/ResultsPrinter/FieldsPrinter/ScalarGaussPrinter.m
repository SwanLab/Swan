classdef ScalarGaussPrinter < ScalarPrinter
    
    properties (Access = private)
        ngaus
        nelem
        gaussDescriptor
    end    
    
    methods (Access = public)
        
        function obj = ScalarGaussPrinter(d)
            obj.ngaus = size(d.fieldValues,1);
            obj.nelem = size(d.fieldValues,3);
            obj.init(d);
            obj.createFieldRepresenter();
            obj.print();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                fieldName = fieldsNames{ifield};
                fieldValue = d.(fieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
        function createFieldRepresenter(obj)
            s.fileID          = obj.fileID;
            s.fieldName       = obj.fieldName;
            s.simulationStr   = obj.simulationStr;
            s.iter            = obj.iter;
            s.fieldType       = obj.fieldType;
            s.fieldPosition   = obj.fieldPosition;
            s.gaussDescriptor = obj.gaussDescriptor;
            s.ngaus           = obj.ngaus;
            s.nelem           = obj.nelem;
            s.nComp           = 1;
            s.fieldValues     = obj.fieldValues;
            obj.fieldRepresenter = GaussFieldPrinter(s);
        end
        
    end
    
end

