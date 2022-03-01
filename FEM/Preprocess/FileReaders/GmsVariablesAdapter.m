classdef GmsVariablesAdapter < handle
    
    properties (Access = private)
       adaptedDataBase 
    end
    
    properties (Access = private)
        nnode
        nfacet
        nelem
        nodeInfo
        connecInfo
        facetsInfo 
    end
    
    methods (Access = public)
       
        function obj = GmsVariablesAdapter(readDataBase)
           obj.init(readDataBase);
        end
        
        function aDB = getAdaptedDB(obj)
            obj.obtainCoordinates();
            obj.obtainConnectivities();
            obj.obtainSetsInfo();
            obj.obtainMasterAndSlaveRelation();
            obj.obtainCorners();
            aDB = obj.adaptedDataBase;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.nnode      = d.nnode;
            obj.nfacet     = d.nfacet;
            obj.nelem      = d.nelem;
            obj.nodeInfo   = d.nodeInfo;
            obj.connecInfo = d.connecInfo;
            obj.facetsInfo = d.facetsInfo;
        end
        
        function obtainCoordinates(obj)
            coordi = obj.nodeInfo;
            obj.adaptedDataBase.coord = zeros(obj.nnode,3);
            obj.adaptedDataBase.coord = coordi(:,1:2);
        end
        
        function obtainSetsInfo(obj)
            connect = obj.connecInfo;
            elementInSet = connect(:,4);
            sets = unique(elementInSet);
            nsets = size(sets,1);
            isElemInThisSet = cell(nsets,1);
            for iset = 1:nsets
                setNum = sets(iset);
                isElemInThisSet{iset} = find(elementInSet == setNum);
            end
            obj.adaptedDataBase.isElemInThisSet = isElemInThisSet;
        end
        
        function obtainConnectivities(obj)
            connect = obj.connecInfo;
            obj.adaptedDataBase.connec = connect(:,1:3);
        end
        
        function obtainMasterAndSlaveRelation(obj)
            coord = obj.adaptedDataBase.coord;
            msRelator = MasterSlaveRelator(coord);
            masterSlave = msRelator.getRelation();
            obj.adaptedDataBase.masterSlave = masterSlave;
        end
        
        function obtainCorners(obj)
            coord = obj.adaptedDataBase.coord;
            cD = CellNodesDescriptor(coord);
            corners = cD.cornerNodes;
            obj.adaptedDataBase.corners = corners;
        end 
        
    end
    
end