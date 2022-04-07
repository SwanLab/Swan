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
            obj.filename = 'CreatedMesh';
            obj.init(cParams);
        end
        
        function computeMeshNodes(obj)
            obj.obtainDimensions();
            obj.computeNodeCoordinates();
            obj.connectNodes();
            obj.obtainMasterSlaveNodes();
            obj.writeFEMreadingArchive();       
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
            obj.div = cParams.divUnit*obj.c;
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
            obj.connec = delaunay(obj.coord);
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
            a = Mesh(s);
            a.plot();
        end
        
        function plotVertices(obj)
            vertexIndex(:,1) = 1:obj.nodes.vert;
            MeshCreator.plotNodes(vertexIndex,obj.coord,'blue')
        end
        
        function plotMasterSlaveNodes(obj)
            masterIndex = obj.masterSlaveIndex(:,1);
            slaveIndex  = obj.masterSlaveIndex(:,2);
            MeshCreator.plotNodes(masterIndex,obj.coord,'green')
            MeshCreator.plotNodes(slaveIndex,obj.coord,'red')
        end
        
    end
    
    methods (Static)
        
        function plotNodes(ind,coord,colorValue)
            b = num2str(ind);
            c = cellstr(b);
            dx = 0.01; dy = 0.01;
            x = coord(ind,1)';
            y = coord(ind,2)';
            t = text(x+dx,y+dy,c);
            set(t,'Color',colorValue)
        end
        
    end

end