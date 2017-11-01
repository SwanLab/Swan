classdef DOF
    %DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem, ?Element, ?Solver}, SetAccess = private)
        ndof
    end
    properties (GetAccess = {?Physical_Problem, ?Element}, SetAccess = private)
        idx
    end
    properties (GetAccess = ?Solver, SetAccess = private)
        vR
        vL
    end
    
    methods (Access = ?Physical_Problem)
        % Constructor
        function obj = DOF(nnode,connec,nunkn,npnod,fixnodes)
            % Compute idx
            lnods = connec';
            for i = 1:nnode
                for j = 1:nunkn
                    obj.idx(nunkn*(i-1)+j,:) = nunkn.*lnods(i,:) - nunkn + j;
                end
            end
            obj.ndof = nunkn*npnod;
            
            % *************************************************************
            if (size(fixnodes,1)>0)
                vR = (fixnodes(:,1)-1)*nunkn + fixnodes(:,2);  % Finds the equation number
                vL = setdiff (1:obj.ndof, vR);
            else
                vL = (1:obj.ndof);
                vR = [];
            end
            
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
            obj.vR = vR;
            obj.vL = vL;
        end
        
    end
    
end

