classdef InputFemFilePrinter < FilePrinter
    
    properties (Access = private)
        connectivities
        coordinates
        isElemInThisSet
        masterSlave
        corners
        nnode
        npnod
        nelem
        ndim
        resultsDir
        scale
        pdim
        ptype
        type
    end
    
    methods (Access = public)
        
        function obj = InputFemFilePrinter(inputData)
            obj.init(inputData);
            obj.createOutPutFolder();
        end
        
        function print(obj)
            obj.openFile();
            obj.printDataProblem();
            obj.printCoordinates();
            obj.printConnectivities();
            obj.printDirichletData();
            obj.printMasterSlave();
            obj.printMaterialSets();
            obj.closeFile();
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,d)
            obj.connectivities  = d.connec;
            obj.coordinates     = d.coord;
            obj.isElemInThisSet = d.isElemInThisSet;
            obj.masterSlave     = d.masterSlave;
            obj.corners         = d.corners;
            obj.scale           = d.scale;
            obj.pdim            = d.pdim;  
            obj.ptype           = d.ptype;
            obj.resultsDir      = d.resultsDir;
            obj.fileName        = d.fileName;     
            obj.type            = d.type;
            obj.nnode           = size(d.connec,2);            
            obj.nelem           = size(d.connec,1);            
            obj.npnod           = size(d.coord,1);            
            obj.ndim            = 2;   
        end
        
        function createOutPutFolder(obj)
            dir = obj.resultsDir;
            if ~exist(dir,'dir')
                mkdir(dir)
                addpath(dir)
            end
        end
        
        function printDataProblem(obj)
            iD = obj.fileID;
            fprintf(iD,'Data_prb = { \n');
            fprintf(iD,['''',obj.type,'''','; \n']);           
            fprintf(iD,'''SI''; \n');                       
            fprintf(iD,['''',obj.pdim,'''','; \n']);           
            fprintf(iD,'''Plane_Stress''; \n');           
            fprintf(iD,['''',obj.ptype,'''','; \n']);           
            fprintf(iD,['''',obj.scale,'''','; \n']); 
            fprintf(iD,'}; \n');            
        end
        
        function printCoordinates(obj)
            iD      = obj.fileID;
            coor    = obj.coordinates;            
            nodInEl = obj.npnod;
            nDim    = obj.ndim;
            fprintf(iD,'\n');
            fprintf(iD,'coord = [\n');
            printFormat = ['%6.0f ',repmat('%12.5d ',1,nDim),' 0 \n'];
            toPrint = [1:nodInEl;coor'];
            fprintf(iD,printFormat,toPrint);
            fprintf(iD,']; \n');
        end        
               
        function printConnectivities(obj)
            iD     = obj.fileID;
            con    = obj.connectivities;
            nNodes = obj.nnode;
            fprintf(iD,'\n');            
            fprintf(iD,'connec = [ \n');
            printFormat = [repmat('%6.0f ',1,nNodes+1),' \n'];
            fprintf(iD,printFormat,[1:obj.nelem;con']);
            fprintf(iD,']; \n');
        end        
        
        function printMasterSlave(obj)
            if ~ isempty(obj.masterSlave)
                iD = obj.fileID;
                fprintf(iD,'\n');
                fprintf(iD,'Master_slave = [\n');
                printFormat = [repmat('%12.0d ',1,2),'\n'];
                fprintf(iD,printFormat,obj.masterSlave');
                fprintf(iD,']; \n');
            end
        end
        
        function printDirichletData(obj)
            iD = obj.fileID;
            fprintf(iD,'\n');            
            fprintf(iD,'dirichlet_data = [\n');
            corn = obj.corners;
            ncorners = size(corn,1);            
            repCorn = repmat(corn,1,2)';
            repDir  = repmat((1:obj.ndim)',ncorners,1);
            ToPrint = [repCorn(:)';repDir(:)'];
            printFormat = [repmat('%12.0f ',1,2),'       0 \n'];
            fprintf(iD,printFormat,ToPrint);
            fprintf(iD,']; \n');  
        end
        
        function printMaterialSets(obj)
            iD = obj.fileID;
            fprintf(iD,'\n');            
            fprintf(iD,'MaterialSets = [\n');
            printFormat = [repmat('%12.0d ',1,2),'\n'];
            for imatSet = 1:numel(obj.isElemInThisSet)
                elemSet = obj.isElemInThisSet{imatSet};
                n       = length(elemSet);
                toPrint = [elemSet';imatSet*ones(1,n)];
                fprintf(iD,printFormat,toPrint);
            end
            fprintf(iD,']; \n');            
        end

    end
    
    
    
end