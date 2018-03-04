classdef Interpolation_variable < Interpolation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = compute(obj,interpolation_geometry,order)
            obj.order = order;
            if strcmp(interpolation_geometry.order,order)==1 % if the interpolations are of the same order copy the properties of the geometry on the variable
                obj.isoparametric = interpolation_geometry.isoparametric;
                obj.xpoints = interpolation_geometry.xpoints;
                obj.T = interpolation_geometry.T;
            else
                obj.compute_xpoints_T(interpolation_geometry,order);
            end
            obj.npnod = length(obj.xpoints(:,1));
        end
        
            function compute_xpoints_T(obj,interp_x,order)
                obj.xpoints = inf*ones(1,3);
                obj.isoparametric = Interpolation_variable.set_isoparametric(interp_x,order);
                nelem = size (interp_x.T);
                inode=1;
                
                
                for inode_variable=1:obj.isoparametric.nnode
                        pos_nodes = num2cell(obj.isoparametric.pos_nodes(inode_variable,1:obj.isoparametric.ndime));
                        shape(inode_variable,:) = cell2mat(interp_x.isoparametric.shape(pos_nodes{:}));
                end
                
                for ielem=1:nelem
                    T_elem=interp_x.T(ielem,:);
                    node_position=1;

                    for inode_variable=1:obj.isoparametric.nnode
                        node=zeros(1,3);
%                         pos_nodes = num2cell(obj.isoparametric.pos_nodes(inode_variable,1:ndime),obj.isoparametric.pos_nodes(inode_variable,2));
%                         shape = interp_x.isoparametric.shape(pos_nodes{:});
                        for inode_mesh=1:interp_x.isoparametric.nnode
                            node= node + shape(inode_variable,inode_mesh)*interp_x.xpoints(T_elem(inode_mesh),:);
                            %             Triangle.shape(Triangle2.pos_nodes(inode_q,1),Triangle2.pos_nodes(inode_q,2))*...
                            
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
        function isoparametric = set_isoparametric(interp_x,order)
            switch interp_x.geometry_type
                
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
                case 'Triangle_Linear_Mass'
                    isoparametric=Triangle_Linear_Mass;
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