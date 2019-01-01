classdef VectorPrinter < FieldPrinter ...
                       & NodalFieldPrinter
    
   
    properties (Access = protected)
        fieldType = 'Vector';
    end      
    
    properties (Access = private)
        fieldComponentName
        ndim = 2;
        field2print
        nodes
        nnode
    end
    
    methods (Access = public)
        
        function obj = VectorPrinter(d)
            obj.init(d);
            obj.computeNodes();
            obj.computeFieldVectorValuesInMatrixForm();
            obj.print();   
        end
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                ifieldName = fieldsNames{ifield};
                fieldValue = d.(ifieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
        function print(obj)
            obj.printResultsLineHeader();
            obj.printComponentNamesLine();
            obj.printValuesLine();
            obj.printFieldLines();
            obj.printEndValuesLine();
        end
        
        function printFieldLines(obj)
            iD = obj.fileID;
            pformat = '%6.0f %12.5d %12.5d \n';
            printV  = [obj.nodes, obj.field2print]';
            fprintf(iD,pformat,printV);
        end
        
        function computeFieldVectorValuesInMatrixForm(obj)
            fV = obj.fieldValues;
            d = obj.ndim;
            fieldV = zeros(obj.nnode,d);
            for idim = 1:d
                dofs = d*(obj.nodes-1)+idim;
                fieldV(:,idim) = fV(dofs);
            end
            obj.field2print = fieldV;
        end
        
        function computeNodes(obj)
            fV = obj.fieldValues;
            d = obj.ndim;
            obj.nnode = round(length(fV)/d);
            obj.nodes(:,1) = 1:obj.nnode;
        end
        
        
        function printComponentNamesLine(obj)
            iD = obj.fileID;
            fC = obj.fieldComponentName;
            d = obj.ndim;
            switch d
                case 2
                    fprintf(iD,'ComponentNames  "%sx", "%sy"\n',fC,fC);
                case 3
                    fprintf(iD,'ComponentNames "%sx", "%sy", "%sz"\n',fC,fC,fC);
            end
        end
        
    end
    
end

