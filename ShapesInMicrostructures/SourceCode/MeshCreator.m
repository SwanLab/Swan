classdef MeshCreator < handle
    
    properties (Access = private)
        c
        theta
        div
        nodes
    end
    
    properties (Access = public)
        filename
        coord
        connec
        masterSlaveIndex
    end
    
    methods (Access = public)
        
        function obj = MeshCreator(cParams)
            obj.init(cParams);
        end
        
        function computeMeshNodes(obj)
            obj.obtainDimensions();
            obj.computeNodeCoordinates();
            obj.connectNodes();
            obj.obtainMasterSlaveNodes();
            %obj.writeFEMreadingArchive();
        end
        
        function drawMesh(obj)
           obj.plotCoordinates();
        end
        
        function plotIndicesOfNodes(obj)
            obj.plotVertices();
            obj.plotMasterSlaveNodes();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.c = cParams.c;
            obj.theta = cParams.theta;
            obj.div = round(cParams.divUnit*obj.c);
            obj.filename = cParams.filename;
        end
        
        function obtainDimensions(obj)
            s.nvert = 2*length(obj.c);
            s.div = obj.div;
            a = NodesCalculator.create(s);
            obj.nodes.vert = a.nvert;
            obj.nodes.bound = a.boundNodes;
            obj.nodes.total = a.totalNodes;
        end
        
        function computeNodeCoordinates(obj)
            s.c = obj.c;
            s.theta = obj.theta;
            s.div = obj.div;
            s.nodes = obj.nodes;
            a = NodeCoordinatesComputer(s);
            a.computeCoordinates();
            obj.coord = a.totalCoord;
        end
        
        function connectNodes(obj)
            s.nodes = obj.nodes;
            s.coord = obj.coord;
            s.theta = obj.theta;
            s.div = obj.div;
            a = NodesConnector(s);
            a.computeConnections();
            obj.connec = a.connec;
        end
        
        function obtainMasterSlaveNodes(obj)
            s.coord = obj.coord;
            s.nodes = obj.nodes;
            s.div = obj.div;
            a = MasterSlaveComputer(s);
            a.computeMasterSlaveNodes();
            obj.masterSlaveIndex = a.masterSlaveIndex;
        end

        function writeFEMreadingArchive(obj)
            s.filename = obj.filename;
            s.coord = obj.coord;
            s.connec = obj.connec;
            s.nodes = obj.nodes;
            s.masterSlaveIndex = obj.masterSlaveIndex;
            a = MatlabFileWriter(s);
            a.write();
        end
     
        function plotCoordinates(obj)
            s.coord = obj.coord;
            s.connec = obj.connec;
            m = Mesh.create(s);
            m.plot();
        end
        
        function plotVertices(obj)
            vertexIndex(:,1) = 1:obj.nodes.vert;
            s.coord = obj.coord;
            s.connec = obj.connec;
            m = Mesh.create(s);
            m.plotNodes(vertexIndex,'blue')
        end
        
        function plotMasterSlaveNodes(obj)
            masterIndex = obj.masterSlaveIndex(:,1);
            slaveIndex  = obj.masterSlaveIndex(:,2);
            s.coord = obj.coord;
            s.connec = obj.connec;
            m = Mesh.create(s);
            m.plotNodes(masterIndex,'green')
            m.plotNodes(slaveIndex,'red')
        end
        
    end

end