classdef GmsReader < FileReader
    
    properties (GetAccess = public, SetAccess = private)
        connec
        coord
        isElemInThisSet
        masterSlave
        corners
    end
    
    properties (Access = private)
        readDataBase
        adaptedDB
    end
        
    methods (Access = public)
        
        function obj = GmsReader(path)
            obj.init(path)
        end
        
        function read(obj)
            obj.openFile();
            obj.readQuantities();
            obj.readNodes();
            obj.readConnec();
            obj.readFacetsConnec();
            obj.closeFile();
            obj.adaptVariables();
        end
        
        function d = getDataBase(obj)
            d = obj.adaptedDB;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,path)
            obj.filePath = path;
        end
        
        function readQuantities(obj)
            var = obj.readLine();
            obj.readDataBase.nnode  = var(1);
            obj.readDataBase.nelem  = var(2);
            obj.readDataBase.nfacet = var(3);
        end
        
        function readNodes(obj)
            for inode=1:obj.readDataBase.nnode
                obj.readDataBase.nodeInfo(inode,:) = obj.readLine();
            end
        end
        
        function readConnec(obj)
            for ielem=1:obj.readDataBase.nelem
                obj.readDataBase.connecInfo(ielem,:) = obj.readLine();
            end
        end
        
        function readFacetsConnec(obj)
            for ielem=1:obj.readDataBase.nfacet
                obj.readDataBase.facetsInfo(ielem,:) = obj.readLine();
            end
        end
        
        function adaptVariables(obj)
            gA = GmsVariablesAdapter(obj.readDataBase);
            obj.adaptedDB = gA.getAdaptedDB();
        end
        
        function var = readLine(obj)
            line = fgetl(obj.fid);
            var = str2num(line);
        end
    end

end