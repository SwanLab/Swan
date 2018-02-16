classdef DOF < handle
    %DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem, ?Element, ?Solver}, SetAccess = private)
        
    end
    properties (GetAccess = {?Physical_Problem, ?Element}, SetAccess = private)
        
    end
    properties (GetAccess = public)
        ndof
        idx
        vC  % Constrainted dofs
        vF  % Free index        
        vD % Diriclet index

    end
    
    methods (Static)
        
        function dof = create(nnode,connec,nunkn,npnod,bc,scale)
            
            switch scale
                case 'MICRO'
                    dof = DOF_micro();
                    [dof.vP,dof.vQ] = dof.compute_periodic_nodes(nunkn,bc);
                case 'MACRO'
                    dof = DOF();
            end
            dof.compute_dofs(nnode,connec,nunkn,npnod,bc)
        end
        
    end
    
    methods

        % Constructor
        function obj = DOF()
            
        end
        
        function compute_dofs(obj,nnode,connec,nunkn,npnod,bc)
            obj.idx = obj.compute_idx(connec,nunkn,nnode);
            obj.ndof = nunkn*npnod;
            obj.vD = obj.compute_diriclet_nodes(nunkn,bc.fixnodes); 
            obj.vC = obj.compute_constrained_nodes();
            obj.vF = setdiff(1:obj.ndof,obj.vC);
        end
        
        function vC = compute_constrained_nodes(obj)
            vC = obj.vD;
        end
        
        
    end
    methods (Static)
        
        function idx = compute_idx(connec,nunkn,nnode)
            idx  = zeros(nnode*nunkn,size(connec,1));
            for i = 1:nnode
                for j = 1:nunkn
                    idx(nunkn*(i-1)+j,:) = nunkn*(connec(:,i) - 1) + j;
                end
            end
            
        end
        
        
        function vD = compute_diriclet_nodes(nunkn,fixnodes)
            if (size(fixnodes,1)>0)
                vD = (fixnodes(:,1)-1)*nunkn + fixnodes(:,2);  % Finds the equation number
            else
                vD = [];
            end
        end
        
        
        
    end
    
    
    
    
end


%     



            % *************************************************************
%             switch linearTriangle.type
%                 case {'TRIANGLE','QUADRILATERAL'}
%                     if (size(fixnodes,1)>0)
%                         vR = (fixnodes(:,1)-1)*nunkn + fixnodes(:,2);  % Finds the equation number
%                         vL = setdiff (1:obj.ndof, vR);
%                     else
%                         vL = (1:obj.ndof);
%                         vR = [];
%                     end
%                 case 'LINEAR_TRIANGLE_MIX'
%                     if (size(fixnodes,1)>0)
%                         vR = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  % Finds the equation number
%                         vL = setdiff (1:obj.ndof, vR);
%                     else
%                         vL = (1:obj.ndof);
%                         vR = [];
%                     end
%                 case 'LINEAR_TRIANGLE_MIX_COUPLED'
%                     vR = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  % Finds the equation number
%                     vL = setdiff (1:obj.ndof, vR);
%                     
%                 case 'HEXAHEDRON'
%                     if (size(fixnodes,1)>0)
%                         vR = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  % Finds the equation number
%                         vL = setdiff (1:obj.ndof, vR);
%                     else
%                         vL = (1:obj.ndof);
%                         vR = [];
%                     end
%                 otherwise
%                     error('No existe es tipo de elemento o no ha sido implementado')
%             end