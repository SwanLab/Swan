classdef GaussHeadPrinter < HeadPrinter
    
    properties (Access = protected)
        gaussDescriptor
        etype
        ngaus
        ndim
        posgp
    end
    
    methods (Access = public)
       
        function print(obj,d)
            obj.init(d);
            obj.printInitialLine();
            obj.printFemMatOoHeader();
            obj.printGaussPointsHeader();            
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
        
        
       function printGaussPointsHeader(obj)
            iD = obj.fileID;
            fprintf(iD,'GaussPoints "%s" Elemtype %s\n',obj.gaussDescriptor,obj.etype);
            fprintf(iD,'Number of Gauss Points: %.0f\n',obj.ngaus);
            fprintf(iD,'Nodes not included\n');
            fprintf(iD,'Natural Coordinates: given\n');
            for igaus = 1:obj.ngaus
                for idime = 1:obj.ndim
                    fprintf(iD,'%12.5d ',obj.posgp(igaus,idime));
                end
                fprintf(iD,'\n');
            end
            fprintf(iD,'End GaussPoints\n');
        end
    end
    
end

