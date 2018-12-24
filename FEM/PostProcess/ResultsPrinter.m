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
    end
    
   
    methods (Access = protected)
        
        function init(obj,fileID,testName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,resultsValues)
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
        end
        
        function print(obj)
            for istep = 1:obj.nsteps
                obj.createFileName(istep);
                obj.openFile();
                obj.printHeader();
                obj.printFemMatOoHeader();
                obj.printGaussPointsHeader();
                obj.printResults(1,istep)
                obj.closeFile();
            end
        end
        
    end
    
    methods (Access = private)
        
        function closeFile(obj)
            fclose(obj.fileID);
        end
        
        function createFileName(obj,iS)
            obj.fileName = fullfile('Output',obj.testName,strcat(obj.fileName,'_','u','_',num2str(iS),'.flavia.msh'));
        end
        
        function openFile(obj)
            obj.fileID = fopen(obj.fileName,'w');            
        end
        
       function printHeader(obj)            
            fprintf(obj.fileID,'GiD Post Results File 1.0\n\n');
       end
        
       function printFemMatOoHeader(obj)
           iD = obj.fileID;
           fprintf(iD,'####################################################\n');
           fprintf(iD,'################# FEM-MAT-OO v.1.0 #################\n');
           fprintf(iD,'####################################################\n');
           fprintf(iD,'\n');
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
    
    methods (Abstract, Access = protected)
       printResults(obj)        
    end
    
end

