classdef Solver_Periodic < Solver
    %AnalyticalSolver Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nunkn
        pnods
    end
    
    
    methods (Access = public)
        
        function obj = Solver_Periodic(obj)
        end
        
        function x = solve(obj,x,LHS,RHS,dof)
            n = dof.ndof;
            nlib = size(obj.pnods(1,:),2);
            pglib = zeros(nlib*obj.nunkn,1);
            qglib = zeros(nlib*obj.nunkn,1);
            for iunkn = 1:obj.nunkn
                index_glib = nlib*(iunkn - 1) + [1:nlib];
                pglib(index_glib,1) = (obj.pnods(1,:)-1)*obj.nunkn + iunkn;
                qglib(index_glib,1) = (obj.pnods(2,:)-1)*obj.nunkn + iunkn;
            end
            
            list = [pglib; qglib; dof.vR];
            free = setdiff([1:1:n],list);
            newLHS = [LHS(free,free) , LHS(free,pglib)+LHS(free,qglib);...
                LHS(pglib,free)+LHS(qglib,free), LHS(pglib,pglib)+LHS(pglib,qglib)+LHS(qglib,pglib)+ LHS(qglib,qglib)];
            
            rhs1 = RHS(pglib)+RHS(qglib);
            rhs = [RHS(free) ; rhs1];
            usol = newLHS \rhs;
            x(free)=usol(1:1:size(free,2));
            x(pglib)=usol(size(free,2)+1:1:size(newLHS,1));
            x(qglib)= x(pglib);
        end
        
        function setSolverVariables(obj,data)
            obj.pnods = data.pnodes;
            obj.nunkn =  data.nunkn;
        end
    end
end
