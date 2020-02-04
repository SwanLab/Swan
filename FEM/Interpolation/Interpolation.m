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
    
    methods (Static, Access = public)
        function interpolation = create(mesh,order)
            switch mesh.geometryType
                case 'LINE'
                    switch order
                        case 'LINEAR'
                            interpolation = Line_Linear(mesh,order);
                        otherwise
                            error('Invalid order for element LINE.');
                    end
                case 'TRIANGLE'
                    switch order
                        case 'LINEAR'
                            interpolation = Triangle_Linear(mesh,order);
                        case 'QUADRATIC'
                            interpolation = Triangle_Quadratic(mesh,order);
                        otherwise
                            error('Invalid order for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'LINEAR'
                            interpolation = Quadrilateral_Bilinear(mesh,order);
                        case 'QUADRATIC'
                            warning('PENDING TO BE TRASFORMED TO INTERPOLATION. SEE TRIANGLE_QUADRATIC AS EXAMPLE')
                            interpolation = Quadrilateral_Serendipity(mesh,order);
                        otherwise
                            error('Invalid order for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    interpolation = Tetrahedra_Linear(mesh,order);
                case 'HEXAHEDRA'
                    interpolation = Hexahedra_Linear(mesh,order);
                otherwise
                    error('Invalid mesh type.')
            end
            
          % if interpolation.nnode ~= size(mesh.connec,2)
           %    interpolation.computeCoordAndConnec(mesh);               
          % end
        end
    end
    
    methods (Access = public)
        
        function obj = Interpolation()

        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,mesh,order)
            switch order 
                case 'LINEAR'
                    obj.xpoints = mesh.coord;
                    obj.T       = mesh.connec;
                otherwise
                    [coord,connec] = obj.computeCoordAndConnec(mesh);
                    obj.xpoints = coord;
                    obj.T       = connec;
            end
            obj.npnod = size(obj.xpoints,1);
            obj.nelem = size(obj.T,1);        
        end
        
    end
    
    methods (Access = private)
        
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
        
%          function [xPoints,Tinterp] = computeCoordAndConnec(obj,mesh)
%             [xPoints,Tinterp] = obj.computeConnectivities(mesh);
% %            obj.T = Tinterp;
% %            obj.xpoints = xPoints;
%             obj.npnod = size(obj.xpoints,1);
%         end
        
        function [xPoints,Tinterp] = computeCoordAndConnec(obj,mesh)
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
        
        function ind = findPointInList(obj,node,xPoints)
            match = true(size(xPoints,1),1);
            for idime = 1:size(node,2)
                match = match & xPoints(:,idime) == node(idime);
            end
            ind = find(match);
        end
        
    end
end
