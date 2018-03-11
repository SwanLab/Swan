classdef Interpolation < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T
        xpoints
        isoparametric
        quadrature
        order
        npnod
        nelem
    end    
    methods
        function obj=Interpolation(mesh)
            nnode = size(mesh.connec,2);
            obj.xpoints = mesh.coord;
            obj.T = mesh.connec;
            obj.npnod = length(obj.xpoints(:,1));
            obj.nelem = length(obj.T(:,1));
            switch mesh.geometryType
                case 'TRIANGLE'
                    switch nnode
                        case 3
                            obj.quadrature = Quadrature_Triangle('LINEAR');
                            obj.isoparametric = Triangle_Linear;
                        case 6
                           obj.quadrature = Quadrature_Triangle('QUADRATIC');
                           obj.isoparametric = Triangle_Quadratic;
                        otherwise
                            error('Invalid nnode for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch nnode
                        case 4
                           obj.quadrature = Quadrature_Quadrilateral('LINEAR');
                           obj.isoparametric = Quadrilateral_Bilinear;
                        case 8
                           obj.quadrature = Quadrature_Quadrilateral('QUADRATIC');
                           obj.isoparametric = Quadrilateral_Serendipity;
                        otherwise
                            error('Invalid nnode for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                   obj.isoparametric = Tetrahedra;
                   obj.quadrature = Quadrature_Tetrahedra('LINEAR');
                case 'HEXAHEDRA'
                   obj.isoparametric = Hexahedra;
                   obj.quadrature = Quadrature_Hexahedra('LINEAR');
                otherwise
                    error('Invalid mesh type.')
            end
        end
        function compute_xpoints_T(obj)
            obj.xpoints = inf*ones(1,3);           
            inode=1;                       
            for inode_variable=1:obj.isoparametric.nnode
                pos_nodes = num2cell(obj.isoparametric.pos_nodes(inode_variable,1:obj.isoparametric.ndime));
                shape(inode_variable,:) = cell2mat(obj.isoparametric.shape(pos_nodes{:}));
            end            
            for ielem=1:obj.nelem
                T_elem=obj.T(ielem,:);
                node_position=1;                
                for inode_variable=1:obj.isoparametric.nnode
                    node=zeros(1,3);
                    for inode_mesh=1:obj.isoparametric.nnode
                        node= node + shape(inode_variable,inode_mesh)*obj.xpoints(T_elem(inode_mesh),:);                        
                    end
                    
                    ind= find(obj.xpoints(:,1)== node(1) & obj.xpoints(:,2)== node(2) & obj.xpoints(:,3)== node(3)); % search if the point is already in the list
                    
                    if isempty(ind)
                        obj.xpoints(inode,:)= node;
                        obj.T(ielem,node_position)= inode;
                        inode=inode+1;
                    else
                        obj.T(ielem,node_position)= ind;
                    end                    
                    node_position=node_position+1;
                end
            end
        end
    end
    methods (Static)
        function isoparametric = set_isoparametric(geometry_type,order)
            switch geometry_type
                
                case 'TRIANGLE'
                    switch order
                        case 'CONSTANT'
                            isoparametric = Triangle_Constant;
                        case 'LINEAR'
                            isoparametric = Triangle_Linear;
                        case 'QUADRATIC'
                            isoparametric = Triangle_Quadratic;
                        otherwise
                            error('Invalid nnode for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'CONSTANT'
                            
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