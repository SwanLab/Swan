classdef VectorNodalPrinter < VectorPrinter ...
        & NodalFieldPrinter
    
    
    properties (Access = private)
        field2print
        nodes
        nnode
    end
    
    methods (Access = public)
        
        function obj = VectorNodalPrinter(d)
            obj.init(d);
            obj.nComp = 2;
            obj.createFormatString();
            obj.computeNodes();
            obj.computeFieldVectorValuesInMatrixForm();
            obj.print();
        end
    end
    
    methods (Access = protected) 
        
        function printFieldLines(obj)
            iD = obj.fileID;
            pformat = '%6.0f %12.5d %12.5d \n';
            printV  = [obj.nodes, obj.field2print]';
            fprintf(iD,pformat,printV);
        end
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                ifieldName = fieldsNames{ifield};
                fieldValue = d.(ifieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
    end
    
    methods (Access = private)
        
        function computeFieldVectorValuesInMatrixForm(obj)
            fV = obj.fieldValues;
            d = obj.nComp;
            fieldV = zeros(obj.nnode,d);
            for idim = 1:d
                dofs = d*(obj.nodes-1)+idim;
                fieldV(:,idim) = fV(dofs);
            end
            obj.field2print = fieldV;
        end
        
        function computeNodes(obj)
            fV = obj.fieldValues;
            d = obj.nComp;
            obj.nnode = round(length(fV)/d);
            obj.nodes(:,1) = 1:obj.nnode;
        end
        
    end
    
end

