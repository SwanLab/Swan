classdef MeshMerger < handle

    properties (Access = private)
        nodesA
        nodesB
        allNodes
        mergedNodes
        mergedMesh
    end
    
    properties (Access = private)
        meshA
        meshB
        isMergedNodeA
        isMergedNodeB
    end
    
    methods (Access = public)
        
        function obj = MeshMerger(cParams)
            obj.init(cParams)
            obj.computeNodesA();
            obj.computeNodesB();
            obj.computeAllNodes();
        end
        
        function m = compute(obj)
            obj.computeMergedNodes();
            obj.computeMergedMesh();
            m = obj.computeCanonicalMergedMesh();
        end
        
        function rNodes = computeRemainingNodes(obj)
            aNodes  = obj.allNodes;
            remainingA = true(obj.meshA.nnodes,1);
            remainingB = ~obj.isMergedNodeB;
            remaining  = [remainingA,remainingB];
            rNodes   = aNodes(remaining);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.meshA = cParams.meshA;
            obj.meshB = cParams.meshB;
            obj.isMergedNodeA = cParams.isMergedNodeA;
            obj.isMergedNodeB = cParams.isMergedNodeB;
        end

        function computeNodesA(obj)
            npnod = obj.meshA.nnodes;
            nodes = (1:npnod)';
            obj.nodesA = nodes;
        end

        function computeNodesB(obj)
            npnodA = obj.meshA.nnodes;
            npnodB = obj.meshB.nnodes;
            nodes = (1:npnodB)' + npnodA;
            obj.nodesB = nodes;
        end

        function computeAllNodes(obj)
            nA = obj.nodesA;
            nB = obj.nodesB;
            nodes = [nA;nB];
            obj.allNodes = nodes;
        end

        function computeMergedNodes(obj)
            nA   = obj.nodesA;
            nB   = obj.nodesB;
            isMa = obj.isMergedNodeA;
            isMb = obj.isMergedNodeB;
            mNodesA          = nA;
            mNodesB(isMb,1)  = nA(isMa);
            mNodesB(~isMb,1) = nB(~isMb);
            nodes = [mNodesA;mNodesB];
            obj.mergedNodes = nodes;
        end

       function computeMergedMesh(obj)
            s.connec = obj.computeMergedConnec();
            s.coord  = obj.computeMergedCoord();
            m = Mesh.create(s);
            obj.mergedMesh = m;
        end

        function mConnec = computeMergedConnec(obj)
            aConnec = obj.computeAllConnec();
            s.oldNodes = obj.allNodes;
            s.newNodes = obj.mergedNodes;
            c = ConnecRenumbering(s);
            mConnec = c.renumber(aConnec);
        end
        
        function connec = computeAllConnec(obj)
            npnodA  = obj.meshA.nnodes;
            connecA = obj.meshA.connec;
            connecB = obj.meshB.connec + npnodA;
            connec = [connecA;connecB];
        end

        function mCoord = computeMergedCoord(obj)
            coordA = obj.meshA.coord;
            coordB = obj.meshB.coord;
            mCoord = [coordA;coordB];
        end

        function m = computeCanonicalMergedMesh(obj)
            s.remainingNodes = obj.computeRemainingNodes();
            s.mesh           = obj.mergedMesh;
            c = CannonicalMeshComputer(s);
            m = c.compute();
        end

    end
    
end