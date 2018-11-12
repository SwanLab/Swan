classdef Interpolation < handle
    properties (GetAccess = public, SetAccess = private)
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
        dvolu
        iteration
        cases
        selectcases
        main_loop
        extra_cases
    end
    
    methods (Static, Access = public)
        function interpolation = create(mesh,order)
            switch mesh.geometryType
                case 'TRIANGLE'
                    switch order
                        case 'LINEAR'
                            interpolation = Triangle_Linear(mesh);
                        case 'QUADRATIC'
                            interpolation = Triangle_Quadratic(mesh);
                        otherwise
                            error('Invalid nnode for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'LINEAR'
                            interpolation = Quadrilateral_Bilinear(mesh);
                        case 'QUADRATIC'
                            warning('PENDING TO BE TRASFORMED TO INTERPOLATION. SEE TRIANGLE_QUADRATIC AS EXAMPLE')
                            interpolation = Quadrilateral_Serendipity(mesh);
                        otherwise
                            error('Invalid nnode for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    interpolation = Tetrahedra_Linear(mesh);
                case 'HEXAHEDRA'
                    interpolation = Hexahedra_Linear(mesh);
                otherwise
                    error('Invalid mesh type.')
            end
            
            if interpolation.nnode ~= size(mesh.connec,2)
                mesh_interpolation = Interpolation.create(mesh,'LINEAR');
                interpolation.compute_xpoints_T(mesh_interpolation)
            end
        end
    end
    
    methods (Access = public)
        function obj = Interpolation(mesh)
            obj.xpoints = mesh.coord;
            obj.T = mesh.connec;
            obj.npnod = size(obj.xpoints,1);
            obj.nelem = size(obj.T,1);
        end
        
        function compute_xpoints_T(obj,mesh_interpolation)
            obj.xpoints = inf*ones(1,size(obj.xpoints,2));
            inode = 1;
            for inode_variable = 1:obj.nnode
                posnodes = obj.pos_nodes(inode_variable,1:obj.ndime);
                mesh_interpolation.computeShapeDeriv(posnodes')
                shape_new(inode_variable,:) = mesh_interpolation.shape;
            end
            
            for ielem = 1:obj.nelem
                T_elem = mesh_interpolation.T(ielem,:);
                node_position = 1;
                for inode_variable = 1:obj.nnode
                    node = zeros(1,size(obj.xpoints,2));
                    for inode_mesh = 1:mesh_interpolation.nnode
                        node =  node + shape_new(inode_variable,inode_mesh)*mesh_interpolation.xpoints(T_elem(inode_mesh),:);
                    end
                    
                    ind = obj.findPointInList(node);
                    
                    if isempty(ind)
                        obj.xpoints(inode,:) =  node;
                        obj.T(ielem,node_position) =  inode;
                        inode = inode+1;
                    else
                        obj.T(ielem,node_position) =  ind;
                    end
                    node_position = node_position+1;
                end
            end
            obj.npnod = size(obj.xpoints,1);
        end
    end
    
    methods (Access = private)
        function ind = findPointInList(obj,node)
            match = true(size(obj.xpoints,1),1);
            for idime = 1:size(node,2)
                match = match & obj.xpoints(:,idime) =  =  node(idime);
            end
            ind = find(match);
        end
    end
end
