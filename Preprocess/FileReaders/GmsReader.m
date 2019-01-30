classdef GmsReader < handle
    
    properties (SetAccess = private)
        connec
        coord
        isElemInThisSet
        masterSlave
        corners
    end 
        
    properties (Access = private)
        filePath
        fid
        nnode
        nfacet
        nelem
        nodeInfo
        connecInfo
        facetsInfo
    end
    
    
    methods (Access = public)
        
        function obj = GmsReader(path)
            obj.init(path)
            obj.readFile();
            obj.obtainCoordinates();
            obj.obtainConnectivities();
            obj.obtainSetsInfo();  
            obj.obtainMasterAndSlaveRelation();
            obj.obtainCorners();
        end
        
    end
    
    methods (Access = private)
    
        function readFile(obj)
            obj.openFile();
            obj.readQuantities();
            obj.readNodes();
            obj.readConnec();
            obj.readFacetsConnec();
            obj.closeFile();
        end
        
        function init(obj,path)
            obj.filePath = path;
        end
        
        function openFile(obj)
            obj.fid = fopen(obj.filePath);
        end
        
        function readQuantities(obj)
            var = obj.readLine();
            obj.nnode  = var(1);
            obj.nelem  = var(2);
            obj.nfacet = var(3);
        end
        
        function readNodes(obj)
            for inode=1:obj.nnode
                obj.nodeInfo(inode,:) = obj.readLine();
            end
        end
        
        function readConnec(obj)
            for ielem=1:obj.nelem
                obj.connecInfo(ielem,:) = obj.readLine();
            end
        end
        
        function readFacetsConnec(obj)
            for ielem=1:obj.nfacet
                obj.facetsInfo(ielem,:) = obj.readLine();
            end
        end
        
        function var = readLine(obj)
            line = fgetl(obj.fid);
            var = str2num(line);
        end
        
        function closeFile(obj)
            fclose(obj.fid);
        end
        
        
        function obtainCoordinates(obj)
            coordi = obj.nodeInfo;
            obj.coord = zeros(obj.nnode,3);
            obj.coord = coordi(:,1:2);
        end
        
        function obtainSetsInfo(obj)
            connect = obj.connecInfo;
            elementInSet = connect(:,4);
            sets = unique(elementInSet);
            for iset = 1:size(sets)
                setNum = sets(iset);
                obj.isElemInThisSet{iset} = find(elementInSet == setNum);
            end
        end
        
        function obtainConnectivities(obj)
            connect = obj.connecInfo;
            obj.connec = connect(:,1:3);
        end
        
        function obtainMasterAndSlaveRelation(obj)
            msRelator = MasterSlaveRelator(obj.coord);
            obj.masterSlave = msRelator.getRelation();
        end
        
        function obtainCorners(obj)
            cD = CellNodesDescriptor(obj.coord);
            obj.corners = cD.cornerNodes;
        end
        
    end
    
    
end