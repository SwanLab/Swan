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
        results
        istep
    end
    
   
    methods (Access = protected)
        
        function init(obj,fileID,testName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,resultsValues,iter)
            obj.fileID   = fileID;
            obj.testName = testName;
            obj.nsteps   = nsteps;
            obj.gaussDescriptor = gaussDescriptor;
            obj.etype    = etype;
            obj.ptype    = ptype;
            obj.ngaus    = ngaus;
            obj.ndim     = ndim;
            obj.posgp    = posgp;
            obj.results  = resultsValues;
            obj.istep    = iter;
        end
        
        function print(obj)
            obj.createFileName();
            obj.openFile();
            obj.printInitialLine();
            obj.printFemMatOoHeader();
            obj.printHeader();
            obj.printResults();
            obj.closeFile();
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

