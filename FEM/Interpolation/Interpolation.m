classdef Interpolation < handle
    
    properties (GetAccess = public, SetAccess = protected)
        T
        xpoints
        ndime
        nnode
        order
        npnod
        nelem
        type
        pos_nodes
        shape
        deriv
        isoDv
        iteration
        cases
        selectcases
        main_loop
        extra_cases
    end
    
    properties (Access = protected)
        mesh
    end
    
    methods (Static, Access = public)
        
        function obj = create(mesh,order)
            cParams.mesh = mesh;
            cParams.order = order;
            f = InterpolationFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.order = cParams.order;            
        end
        
        function computeCoordAndConnec(obj)
            switch obj.order
                case 'LINEAR'
                    coord  = obj.mesh.coord;
                    connec = obj.mesh.connec;
                otherwise
                    [coord,connec] = obj.computeCoordAndConnecInGeneral(obj.mesh);                    
            end
            obj.xpoints = coord;
            obj.T       = connec;
            obj.npnod = size(obj.xpoints,1);
            obj.nelem = size(obj.T,1);
        end
        
    end
    
    methods (Access = private)
        
        function [xPoints,Tinterp] = computeCoordAndConnecInGeneral(obj,mesh)
            shapesInVarNodes = obj.computeShapesInVariableNodes(mesh);
            xPointsMesh      = mesh.coord;
            Tmesh            = mesh.connec;
            
            xPoints = zeros(1,obj.ndime);
            Tinterp = zeros(obj.nelem,obj.nnode);
            
            inode = 1;
            for ielem = 1:mesh.nelem
                for inodeVar = 1:obj.nnode
                    xNode = zeros(1,obj.ndime);
                    for inodeMesh = 1:mesh.nnode
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
        end
        
        function shapes = computeShapesInVariableNodes(obj,mesh)
            interpMesh = Interpolation.create(mesh,'LINEAR');
            nNodeMesh = interpMesh.nnode;
            nNodeVar  = obj.nnode;
            shapes    = zeros(nNodeVar,nNodeMesh);
            nodesVar  = obj.pos_nodes;
            for inodeVar = 1:obj.nnode
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
