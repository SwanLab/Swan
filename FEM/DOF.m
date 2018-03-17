classdef DOF < handle
    %DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem, ?Element, ?Solver}, SetAccess = private)
        
    end
    properties (GetAccess = {?Physical_Problem, ?Element}, SetAccess = private)
        
    end
    properties (GetAccess = public)
        dirichlet_values
        neumann_values
        full_dirichlet_values
        ndof
        nunkn
        in_elem
        constrained  % Constrainted dofs index
        free  % Free dof index
        dirichlet % Diriclet dof index
        neumann % Explicit (from input) Neumann dof
        full_dirichlet % Everywhere dirichlet dof inex
        periodic_free % Perioic
        periodic_constrained
    end
    
    methods
        function obj = DOF(filename,geometry,mesh)
            switch mesh.ptype                
                case 'ELASTIC'
                    switch mesh.pdim
                        case '2D'
                            obj.nunkn = 2;
                        case '3D'                            
                            obj.nunkn = 3;
                    end
                case 'THERMAL'
                    obj.nunkn=1;
                case 'DIFF-REACT'
                    obj.nunkn=1;                    
                case 'HYPERELASTIC'
                    obj.nunkn=2;                 
            end    
            switch mesh.ptype
                case 'Stokes'
                    nunkn_u = 2;
                    nunkn_p = 1;
                    obj.nunkn = [nunkn_u nunkn_p];
                    [dirichlet_data,neumann_data,full_dirichlet_data,master_slave] = ...
                        Preprocess.getBC_fluids(filename,geometry);
                otherwise
                    [dirichlet_data,neumann_data,full_dirichlet_data,master_slave] = Preprocess.getBC_mechanics(filename);
            end
            
            
            [obj.neumann,obj.neumann_values] = obj.get_dof_conditions(neumann_data,obj.nunkn(1));
            [obj.full_dirichlet,obj.full_dirichlet_values] = obj.get_dof_conditions(full_dirichlet_data,obj.nunkn(1));
            
           for ifield = 1:geometry(1).nfields
                nunkn = obj.nunkn(ifield);
                nnode = geometry(ifield).interpolation.isoparametric.nnode;
                npnod = geometry(ifield).interpolation.npnod;
                obj.in_elem{ifield} = obj.compute_idx(geometry(ifield).interpolation.T,nunkn,nnode);


                obj.ndof(ifield) = nunkn*npnod;

                [obj.dirichlet{ifield},obj.dirichlet_values{ifield}] = obj.get_dof_conditions(dirichlet_data{ifield},nunkn);
                
                

                if ~isempty(master_slave)
                obj.periodic_free = obj.compute_periodic_nodes(master_slave(:,1),nunkn);
                obj.periodic_constrained = obj.compute_periodic_nodes(master_slave(:,2),nunkn);
                end

                obj.constrained{ifield} = obj.compute_constrained_dof(mesh.scale,ifield);
                obj.free{ifield} = obj.compute_free_dof(ifield);
            end
%                 if nfields == 1
%                     obj.dirichlet = cell2mat(obj.dirichlet{ifield});
%                     obj.dirichlet_values = cell2mat(obj.dirichlet_values{ifield});
%                     
%                 end
        end
        
    end
    
    methods
        
        function constrained = compute_constrained_dof(obj,scale,ifield)
            switch scale
                case 'MICRO'
                    constrained = [obj.periodic_constrained;obj.dirichlet{ifield}];
                case 'MACRO'
                    constrained = obj.dirichlet{ifield};
            end
        end
        
        function free = compute_free_dof(obj,ifield)
            free = setdiff(1:obj.ndof(ifield),obj.constrained{ifield});
        end
        
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
        
        
        function periodic_dof = compute_periodic_nodes(obj,periodic_nodes,nunkn)
            nlib = size(periodic_nodes,1);
            periodic_dof = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                index_glib = nlib*(iunkn - 1) + [1:nlib];
                periodic_dof(index_glib,1) = obj.unknown_and_node_id_to_dof_id(periodic_nodes,iunkn,nunkn);
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
        
        function  [dof_condition_id,dof_condition_value] = get_dof_conditions(obj,conditions_unkn_and_dim,nunkn)
            if isempty(conditions_unkn_and_dim)
                dof_condition_id = [];
                dof_condition_value = [];
            else
          
                if ismatrix(conditions_unkn_and_dim)
                    inode_condition = conditions_unkn_and_dim(:,1);
                    iunkn_condition = conditions_unkn_and_dim(:,2);
                    dof_condition_id = obj.unknown_and_node_id_to_dof_id(inode_condition,iunkn_condition,nunkn);
                    dof_condition_value = conditions_unkn_and_dim(:,3);
                else
                    dof_condition_id = [];
                    dof_condition_value = conditions_unkn_and_dim;
                end
            end
        end
        
        
    end
    
    methods (Static)
        
        function idof = unknown_and_node_id_to_dof_id(inode,iunkn,nunkn)
            idof(:,1)= nunkn*(inode - 1) + iunkn;
        end
        
        
        function [Master_slave_nodes] = get_master_slave_in_square(coordinates)
            % Square muest be [0,1]x[0,1]
            % nodes in the left-vertical side, without the corners
            href = 0.025;
            h=href; L=[];Y=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,1)<h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5 )
                    L = [L i];
                    Y = [Y coordinates(i,2)];
                end
            end
            [Y1,I] = sort(Y);
            V1 = L(I);
            
            % nodes in the right-vertical side, without the corners
            h=href; L=[];Y=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,1)>1-h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5)
                    L = [L i];
                    Y = [Y coordinates(i,2)];
                end
            end
            [Y1,I] = sort(Y);
            V2 = L(I);
            
            % nodes in the bottom-horizontal side, without the corners
            h=href; L=[];X=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,2)<h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5 )
                    L = [L i];
                    X = [X coordinates(i,1)];
                end
            end
            [X1,I] = sort(X);
            H1 = L(I);
            
            % nodes in the top-horizontal side, without the corners
            h=href; L=[];X=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,2)>1-h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5)
                    L = [L i];
                    X = [X coordinates(i,1)];
                end
            end
            [X1,I] = sort(X);
            H2 = L(I);
            Master_slave_nodes = [V1 H1; V2 H2]; % lista de nodos
        end
        
        
    end
    
    
end

