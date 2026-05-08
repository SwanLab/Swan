% classdef MeshCreator < handle
% 
%     properties (Access = private)
%         c
%         theta
%         div
%         nodes
%         latticeVectors
%     end
% 
%     properties (Access = public)
%         filename
%         coord
%         connec
%         masterSlaveIndex
%         lattice
%     end
% 
%     methods (Access = public)
% 
%         function obj = MeshCreator(cParams)
%             obj.init(cParams);
%         end
% 
%         function computeMeshNodes(obj)
%             if ~isempty(obj.latticeVectors)
%                 obj.computeCAndThetaFromVectors();
%             end
%             obj.obtainDimensions();
%             obj.computeNodeCoordinates();
%             obj.connectNodes();
%             obj.obtainMasterSlaveNodes();
% %             obj.writeFEMreadingArchive();
% 
%             % v = obj.coord(1:obj.nodes.vert,:);
%             % a1 = v(2,:) - v(1,:);
%             % a2 = v(3,:) - v(2,:);
%             % obj.lattice = LatticeVectors(a1,a2);
%         end
% 
%         function drawMesh(obj)
%            obj.plotCoordinates();
%         end
% 
%         function plotIndicesOfNodes(obj)
%             obj.plotVertices();
%             obj.plotMasterSlaveNodes();
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function init(obj,cParams)
%             if isfield(cParams,'latticeVectors')
%                 obj.latticeVectors = cParams.latticeVectors;
%                 obj.c = [];
%                 obj.theta = [];
%             else
%                 obj.c = cParams.c;
%                 obj.theta = cParams.theta;
%                 obj.latticeVectors = [];
%             end
%             if ~isempty(obj.c)
%                 obj.div = round(cParams.divUnit * obj.c);
%             else
%                 lengths = vecnorm(obj.latticeVectors, 2, 2)';
%                 obj.div = round(cParams.divUnit * lengths);
%             end
%             obj.filename = cParams.filename;
%         end
% 
%          function computeCAndThetaFromVectors(obj)
%              nVectors = size(obj.latticeVectors, 1);
%              obj.c = zeros(1, nVectors);
%              obj.theta = zeros(1, nVectors);
%              for i = 1:nVectors
%                  v = obj.latticeVectors(i, :);
%                  obj.c(i) = norm(v);
%                  obj.theta(i) = atan2(v(2), v(1)) * 180/pi;
%              end
%          end
% 
%         function obtainDimensions(obj)
%             s.nvert = 2*length(obj.c);
%             s.div = obj.div;
%             a = NodesCalculator.create(s);
%             obj.nodes.vert = a.nvert;
%             obj.nodes.bound = a.boundNodes;
%             obj.nodes.total = a.totalNodes;
%         end
% 
%         function computeNodeCoordinates(obj)
%             s.c = obj.c;
%             s.theta = obj.theta;
%             s.div = obj.div;
%             s.nodes = obj.nodes;
%             a = NodeCoordinatesComputer(s);
%             a.computeCoordinates();
%             obj.coord = a.totalCoord;
%         end
% 
%         function connectNodes(obj)
%             s.nodes = obj.nodes;
%             s.coord = obj.coord;
%             s.theta = obj.theta;
%             s.div = obj.div;
%             a = NodesConnector(s);
%             a.computeConnections();
%             obj.connec = a.connec;
%         end
% 
%         function obtainMasterSlaveNodes(obj)
%             s.coord = obj.coord;
%             s.nodes = obj.nodes;
%             s.div = obj.div;
%             a = MasterSlaveComputer(s);
%             a.computeMasterSlaveNodes();
%             obj.masterSlaveIndex = a.masterSlaveIndex;
%         end
% 
%         function writeFEMreadingArchive(obj)
%             s.filename = obj.filename;
%             s.coord = obj.coord;
%             s.connec = obj.connec;
%             s.nodes = obj.nodes;
%             s.masterSlaveIndex = obj.masterSlaveIndex;
%             a = MatlabFileWriter(s);
%             a.write();
%         end
% 
%         function plotCoordinates(obj)
%             s.coord = obj.coord;
%             s.connec = obj.connec;
%             m = Mesh.create(s);
%             m.plot();
%         end
% 
%         function plotVertices(obj)
%             vertexIndex(:,1) = 1:obj.nodes.vert;
%             s.coord = obj.coord;
%             s.connec = obj.connec;
%             m = Mesh.create(s);
%             m.plotNodes(vertexIndex,'blue')
%         end
% 
%         function plotMasterSlaveNodes(obj)
%             masterIndex = obj.masterSlaveIndex(:,1);
%             slaveIndex  = obj.masterSlaveIndex(:,2);
%             s.coord = obj.coord;
%             s.connec = obj.connec;
%             m = Mesh.create(s);
%             m.plotNodes(masterIndex,'green')
%             m.plotNodes(slaveIndex,'red')
%         end
% 
%     end
% 
% end

classdef MeshCreator < handle
    
    properties (Access = private)
        latticeVectors
        div
        nodes
    end
    
    properties (Access = public)
        filename
        coord
        connec
        masterSlaveIndex
        lattice
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
            obj.latticeVectors = cParams.latticeVectors;            
            
            nVectors = size(obj.latticeVectors, 1);
            obj.div = repmat(cParams.divUnit, 1, nVectors);
            
            obj.filename = cParams.filename;
        end
        
        function obtainDimensions(obj)
            nVectors = size(obj.latticeVectors, 1);
            s.nvert = 2 * nVectors;
            s.div = obj.div;
            a = NodesCalculator.create(s);
            obj.nodes.vert = a.nvert;
            obj.nodes.bound = a.boundNodes;
            obj.nodes.total = a.totalNodes;
        end
        
        function computeNodeCoordinates(obj)
            s.latticeVectors = obj.latticeVectors;
            s.div = obj.div;
            s.nodes = obj.nodes;
            a = NodeCoordinatesComputer(s);
            a.computeCoordinates();
            obj.coord = a.totalCoord;
        end
        
        function connectNodes(obj)
            s.nodes = obj.nodes;
            s.coord = obj.coord;
            s.latticeVectors = obj.latticeVectors;
            s.div = obj.div;
            a = NodesConnector(s);
            a.computeConnections();
            obj.connec = a.connec;
        end
        
        function obtainMasterSlaveNodes(obj)
            s.coord = obj.coord;
            s.nodes = obj.nodes;
            s.div = obj.div;
            s.latticeVectors = obj.latticeVectors;
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