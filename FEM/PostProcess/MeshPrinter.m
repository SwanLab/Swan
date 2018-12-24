classdef MeshPrinter < handle
    
    properties (Access = private)
        nsteps
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
    end
    
    methods (Access = public)
        
        function obj = MeshPrinter(nsteps,testName,npnod,pdim,nnode,coordinates,connectivities,nelem,ndim,etype)
            obj.init(nsteps,testName,npnod,pdim,nnode,coordinates,connectivities,nelem,ndim,etype)
            obj.print()
        end
    end
    
    methods (Access = private)
        
        
        function init(obj,nsteps,testName,npnod,pdim,nnode,coordinates,connectivities,nelem,ndim,etype)
            obj.nsteps    = nsteps;
            obj.testName  = testName;
            obj.npnod     = npnod;
            obj.pdim      = pdim;
            obj.nnode     = nnode;
            obj.coordinates = coordinates;
            obj.connectivities = connectivities;
            obj.nelem = nelem;
            obj.ndim = ndim;
            obj.etype = etype;
        end
        
        function print(obj)
            for istep = 1:obj.nsteps
                obj.createFileName(istep);
                obj.openFile();
                obj.printFemMatOoHeader();
                obj.printHeader();
                obj.printCoordinates();
                obj.printConnectivities();
                obj.closeFile();
            end
        end
        
        function createFileName(obj,iS)
            obj.fileName = fullfile('Output',obj.testName,strcat(obj.fileName,'_','u','_',num2str(iS),'.flavia.msh'));
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

