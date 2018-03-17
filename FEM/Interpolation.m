classdef Interpolation < handle
    %Interpolation Summary of this class goes here
    % Detailed explanation goes here
    
    properties
        T
        xpoints
        isoparametric
        order
        npnod
        nelem
    end
    methods
        function obj = Interpolation(mesh,order)
            nnode = size(mesh.connec,2);
            obj.xpoints = mesh.coord;
            obj.T = mesh.connec;
            obj.npnod = length(obj.xpoints(:,1));
            obj.nelem = length(obj.T(:,1));
            obj.isoparametric = obj.createIsoparametric(mesh.geometryType,order);
            if nnode ~= obj.isoparametric.nnode
                mesh_isoparametric = obj.createIsoparametric(mesh.geometryType,'LINEAR');
                obj.compute_xpoints_T(nnode,obj.xpoints,obj.T,mesh_isoparametric.shape)
            end
        end
        function compute_xpoints_T(obj,nnode,xpoints,connec,shape_ini)
            obj.xpoints = inf*ones(1,3);
            inode = 1;
            for inode_variable = 1:obj.isoparametric.nnode
                pos_nodes = num2cell(obj.isoparametric.pos_nodes(inode_variable,1:obj.isoparametric.ndime));
                shape(inode_variable,:) = cell2mat(shape_ini(pos_nodes{:}));
            end
            for ielem = 1:obj.nelem
                T_elem = connec(ielem,:);
                node_position = 1;
                for inode_variable = 1:obj.isoparametric.nnode
                    node = zeros(1,3);
                    for inode_mesh = 1:nnode
                        node = node + shape(inode_variable,inode_mesh)*xpoints(T_elem(inode_mesh),:);
                    end
                    
                    ind = find(obj.xpoints(:,1) == node(1) & obj.xpoints(:,2) == node(2) & obj.xpoints(:,3) == node(3)); % search if the point is already in the list
                    
                    if isempty(ind)
                        obj.xpoints(inode,:) = node;
                        obj.T(ielem,node_position) = inode;
                        inode = inode+1;
                    else
                        obj.T(ielem,node_position) = ind;
                    end
                    node_position = node_position+1;
                end
            end
            obj.npnod = length(obj.xpoints(:,1));
        end
    end
    methods (Static)
        function isoparametric = createIsoparametric(type,order)
            switch type
                case 'TRIANGLE'
                    switch order
                        case 'LINEAR'
                            isoparametric = Triangle_Linear;
                        case 'QUADRATIC'
                            isoparametric = Triangle_Quadratic;
                        otherwise
                            error('Invalid nnode for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'LINEAR'
                            isoparametric = Quadrilateral_Bilinear;
                        case 'QUADRATIC'
                            isoparametric = Quadrilateral_Serendipity;
                        otherwise
                            error('Invalid nnode for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    isoparametric = Tetrahedra;
                case 'HEXAHEDRA'
                    isoparametric = Hexahedra;
                otherwise
                    error('Invalid mesh type.')
            end
        end
    end
end