classdef ConnecCoordFromInterpAndMesh < handle
    
    properties (Access = public)
        coord
        connec
    end
    
    properties (Access = private)
        mesh
        interp
    end
    
    methods (Access = public)
        
        function obj = ConnecCoordFromInterpAndMesh(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            m = obj.mesh;
            shapesInVarNodes = obj.computeShapesInVariableNodes(m);
            xPointsMesh      = m.coord;
            Tmesh            = m.connec;
            
            xPoints = zeros(1,obj.interp.ndime);
            Tinterp = zeros(m.nelem,obj.interp.nnode);
            
            inode = 1;
            for ielem = 1:m.nelem
                for inodeVar = 1:obj.interp.nnode
                    xNode = zeros(1,obj.interp.ndime);
                    for inodeMesh = 1:m.nnode
                        node = Tmesh(ielem,inodeMesh);
                        shapes = shapesInVarNodes(inodeVar,inodeMesh);
                        xNode = xNode + shapes*xPointsMesh(node,:);
                    end
                    
                    node = obj.findPointInList(xNode,xPoints);
                    
                    if isempty(node)
                        xPoints(inode,:) = xNode;
                        Tinterp(ielem,inodeVar) = inode;
                        inode = inode+1;
                    else
                        Tinterp(ielem,inodeVar) = node;
                    end
                end
            end
            obj.coord = xPoints;
            obj.connec = Tinterp;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.interp = cParams.interpolation;
        end
        
        
        function shapes = computeShapesInVariableNodes(obj,mesh)
            interpMesh = Interpolation.create(mesh,'LINEAR');
            nNodeMesh = interpMesh.nnode;
            nNodeVar  = obj.interp.nnode;
            shapes    = zeros(nNodeVar,nNodeMesh);
            nodesVar  = obj.interp.pos_nodes;
            for inodeVar = 1:obj.interp.nnode
                nodesPoints = nodesVar(inodeVar,:);
                interpMesh.computeShapeDeriv(nodesPoints')
                shapes(inodeVar,:) = interpMesh.shape;
            end
        end
        
        
    end
    
    methods (Access = private, Static)
        
        function ind = findPointInList(node,xPoints)
            match = true(size(xPoints,1),1);
            for idime = 1:size(node,2)
                match = match & xPoints(:,idime) == node(idime);
            end
            ind = find(match);
        end
        
    end
    
    
    
    
    
    
end
