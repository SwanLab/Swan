classdef GaussFieldPrinter < FieldRepresenterPrinter

    properties (Access = private)
        gaussDescriptor
        ngaus
        nelem
    end
    
    methods (Access = public)
        
        function obj = GaussFieldPrinter(cParams)
            obj.init(cParams)
        end
        
        function printResultsLineHeader(obj)
            iD = obj.fileID;
            fN = obj.fieldName;
            sC = obj.simulationStr;
            is = obj.iter;
            fT = obj.fieldType;
            rL = obj.fieldPosition;
            gD = obj.gaussDescriptor;
            fprintf(iD,'\nResult "%s" "%s" %.0f %s %s "%s"\n',fN,sC,is,fT,rL,gD);
        end
        
        function printFieldLines(obj)
            iD = obj.fileID;
            pformat = ['%s ',repmat('%12.5d ',1,obj.nComp),'\n'];
            elemColum = obj.computeElementStringColum();
            valColums = obj.computeTensorValueColums();
            str = [elemColum,valColums]';
            fprintf(iD,pformat,str{:});
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
        
        function c = computeElementStringColum(obj)
            allElem(:,1) = 1:obj.nelem;
            colWidth = size(num2str(obj.nelem),2);
            strInCol = repmat(' ',obj.nelem*obj.ngaus,colWidth);
            numIndex = 1:obj.ngaus:obj.nelem*obj.ngaus;
            strInCol(numIndex,:) = num2str(allElem);
            c = cellstr(strInCol);
        end
        
        function fM = computeTensorValueColums(obj)
            fV = obj.fieldValues;
            fM  = zeros(obj.ngaus*obj.nelem,obj.nComp);
            for istre = 1:obj.nComp
                for igaus = 1:obj.ngaus
                    rows = linspace(igaus,(obj.nelem - 1)*obj.ngaus + igaus,obj.nelem);
                    fM(rows,istre) = fV(igaus,istre,:);
                end
            end
            fM = num2cell(fM);
        end
        
    end
    
end

