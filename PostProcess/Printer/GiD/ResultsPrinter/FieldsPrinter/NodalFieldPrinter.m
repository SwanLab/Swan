classdef NodalFieldPrinter < FieldRepresenterPrinter
    
    properties (Access = private)
        nodes
    end
    
    methods (Access = public)
        
        function obj = NodalFieldPrinter(cParams)
            obj.init(cParams);
            obj.nodes(:,1) = 1:size(obj.fieldValues,1);
        end
        
        function printResultsLineHeader(obj)
            iD = obj.fileID;
            fN = obj.fieldName;
            sC = obj.simulationStr;
            it = obj.iter;
            fT = obj.fieldType;
            rL = obj.fieldPosition;
            fprintf(iD,'\nResult "%s" "%s" %.0f %s %s\n',fN,sC,it,fT,rL);
        end
        
        function printFieldLines(obj)
            iD = obj.fileID;
            pformat = ['%6.0f ',repmat('%12.5d ',1,obj.nComp),' \n'];
            printV  = [obj.nodes, obj.fieldValues]';
            fprintf(iD,pformat,printV);
        end
        
    end
    
     methods (Access = private)
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                ifieldName = fieldsNames{ifield};
                fieldValue = d.(ifieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
     end
    
end

