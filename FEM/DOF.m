classdef DOF < handle
    %DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem, ?Element, ?Solver}, SetAccess = private)
        
    end
    properties (GetAccess = {?Physical_Problem, ?Element}, SetAccess = private)
        
    end
    properties (GetAccess = public)
        dirichlet_nodes
        boundary_nodes
        neumann_nodes
        periodic_nodes
        ndof
        in_elem
        constrained  % Constrainted dofs index
        free  % Free dof index
        dirichlet % Diriclet dof index
        boundary % Boundary dof inex
        neumann % Explicit (from input) Neumann dof
        periodic_free % Perioic
        periodic_constrained
    end
    
    methods
        
        function obj = DOF(filename,nnode,connec,nunkn,npnod,scale)
            
            obj.in_elem = obj.compute_idx(connec,nunkn,nnode);
            [obj.dirichlet_nodes,obj.boundary_nodes,obj.neumann_nodes] = Preprocess.getBC(filename);
            
            obj.ndof = nunkn*npnod;
            
            % Dirichlet
            inode_dirichlet = obj.dirichlet_nodes(:,1);
            iunkn_dirichlet = obj.dirichlet_nodes(:,2);
            obj.dirichlet = obj.unknown_and_node_id_to_dof_id(inode_dirichlet,iunkn_dirichlet,nunkn);
            
            % Neumann
            inode_neumann = obj.neumann_nodes(:,1);
            iunkn_neumann = obj.neumann_nodes(:,2);
            obj.neumann = obj.unknown_and_node_id_to_dof_id(inode_neumann,iunkn_neumann,nunkn);
            
            switch scale
                case 'MICRO'
                    [obj.periodic_free,obj.periodic_constrained] = obj.compute_periodic_nodes(nunkn);
                    obj.periodic_nodes = Preprocess.getPeriodicBC(coords);
                    obj.constrained = [obj.periodic_constrained;obj.dirichlet];
                case 'MACRO'
                    obj.constrained = obj.dirichlet;
            end
            obj.free = setdiff(1:obj.ndof,obj.constrained);
        end
        
    end
    
    methods
        
        % Constructor
        
        
        %         function obj=computeFixedNodesValues(obj,ptype,ndim)
        %             switch ptype
        %                 case 'ELASTIC'
        %                     ifix=1;
        %                     for i=1:size(obj.fixnodes_perimeter,1)
        %                         for j=1:ndim
        %                             obj.fixnodes(ifix,1)=obj.fixnodes_perimeter(i,1); % node
        %                             obj.fixnodes(ifix,2)=j; % idim
        %                             obj.fixnodes(ifix,3)=0; % U_imp
        %                             ifix=ifix+1;
        %                         end
        %                     end
        %             end
        %         end
        
        
        function [periodic_free,periodic_constrained] = compute_periodic_nodes(obj,nunkn)
            nlib = size(obj.periodic_nodes(1,:),2);
            periodic_free = zeros(nlib*nunkn,1);
            periodic_constrained = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                index_glib = nlib*(iunkn - 1) + [1:nlib];
                periodic_free(index_glib,1) = obj.unknown_and_node_id_to_dof_id(obj.periodic_nodes(1,:),iunkn,nunkn);
                periodic_constrained(index_glib,1) = obj.unknown_and_node_id_to_dof_id(obj.periodic_nodes(2,:),iunkn,nunkn);

            end
        end
        
        function dof_elem = compute_idx(obj,connec,nunkn,nnode)
            dof_elem  = zeros(nnode*nunkn,size(connec,1));
            for inode = 1:nnode
                for iunkn = 1:nunkn
                    idof_elem = obj.unknown_and_node_id_to_dof_id(inode,iunkn,nunkn);
                    global_node = connec(:,inode);
                    idof_global = obj.unknown_and_node_id_to_dof_id(global_node,iunkn,nunkn);
                    dof_elem(idof_elem,:) = idof_global;
                end
            end
            
        end
        
        
    end
    
    methods (Static)
        
        function idof = unknown_and_node_id_to_dof_id(inode,iunkn,nunkn)
            idof(:,1)= nunkn*(inode - 1) + iunkn;
        end
    end
    
end

