classdef ResultsPrinter < handle
    
    properties
        fileID
        testName
        fileName
        nsteps
        gaussDescriptor
        etype
        ptype
        ngaus
        ndim
        posgp
        fields
        istep
    end
    
    methods (Access = public)

        function print(obj,d)
            obj.init(d)
            obj.createFileName();
            obj.openFile();
            obj.printInitialLine();
            obj.printFemMatOoHeader();
            obj.printHeader();
            obj.printResults();
            obj.closeFile();
        end 
        
        function f = getFieldName(obj)
            f = obj.fileName;
        end
        
    end
   
    methods (Access = protected)
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                fieldName = fieldsNames{ifield};
                fieldValue = d.(fieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
    end
    
    methods (Access = private)
        
        function closeFile(obj)
            fclose(obj.fileID);
        end
        
        function createFileName(obj)
            iS = obj.istep;
            obj.fileName = fullfile('Output',obj.testName,strcat(obj.testName,num2str(iS),'.flavia.res'));
        end
        
        function openFile(obj)
            obj.fileID = fopen(obj.fileName,'w');            
        end
        
       function printInitialLine(obj)            
            fprintf(obj.fileID,'GiD Post Results File 1.0\n\n');
       end
       
        
       function printFemMatOoHeader(obj)
           iD = obj.fileID;
           fprintf(iD,'####################################################\n');
           fprintf(iD,'################# FEM-MAT-OO v.1.0 #################\n');
           fprintf(iD,'####################################################\n');
           fprintf(iD,'\n');
       end
               
    end
    
    methods (Abstract, Access = protected)
       printResults(obj) 
       printHeader(obj)
    end
    
end

