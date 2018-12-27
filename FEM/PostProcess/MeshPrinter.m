classdef MeshPrinter < handle
    
    properties (Access = private)
        testName
        fileID
        npnod
        pdim
        nnode
        coordinates
        connectivities
        nelem
        ndim
        etype
        fileName
        iter
    end
    
    methods (Access = public)
        
        function obj = MeshPrinter()
        end
        
        function print(obj,d)
            obj.init(d);
            obj.createFileName();
            obj.openFile();
            obj.printFemMatOoHeader();
            obj.printHeader();
            obj.printCoordinates();
            obj.printConnectivities();
            obj.closeFile();
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
              
        function createFileName(obj)
            iS = obj.iter;
            obj.fileName = fullfile('Output',obj.testName,strcat(obj.testName,num2str(iS),'.flavia.msh'));
        end
        
        function openFile(obj)
            obj.fileID = fopen(obj.fileName,'w');
        end
        
        function printFemMatOoHeader(obj)
            iD = obj.fileID;
            fprintf(iD,'####################################################\n');
            fprintf(iD,'################# FEM-MAT-OO v.1.0 #################\n');
            fprintf(iD,'####################################################\n');
            fprintf(iD,'\n');
        end
        
        function printCoordinates(obj)
            iD     = obj.fileID;
            coord  = obj.coordinates;
            nodInEl = obj.npnod;
            nDim   = obj.ndim;
            fprintf(iD,'coordinates \n');
            printFormat = ['%6.0f ',repmat('%12.5d ',1,nDim),'\n'];
            toPrint = [1:nodInEl;coord'];
            fprintf(iD,printFormat,toPrint);
            fprintf(iD,'end coordinates \n \n');
        end
        
        function printConnectivities(obj)
            iD     = obj.fileID;
            connec = obj.connectivities;
            nNodes = obj.nnode;
            fprintf(iD,'elements \n');
            printFormat = [repmat('%6.0f ',1,nNodes+1),' 1 \n'];
            fprintf(iD,printFormat,[1:obj.nelem;connec']);
            fprintf(iD,'end elements\n\n');
        end
        
        function printHeader(obj)
            iD = obj.fileID;
            nD = obj.ndim;
            eT = obj.etype;
            nN = obj.nnode;
            printFormat = 'MESH "WORKPIECE" dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n';
            fprintf(iD,printFormat,nD,eT,nN);
        end
        
        function closeFile(obj)
            fclose(obj.fileID);
        end
        
    end
    
    
    
end

