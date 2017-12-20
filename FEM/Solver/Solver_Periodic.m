classdef Solver_Periodic < Solver
    %AnalyticalSolver Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Access = {?Physical_Problem,?Physical_Problem_Micro})
        function x = solve(LHS,RHS,dof,fixnodes,pnods)
            x = zeros(dof.ndof,1);
            n = size(dof.ndof,1);
            nlib = size(pnods(1,:),2);%%%%%%
            pglib = zeros(nlib*nunkn,1);
            qglib = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                index_glib = nlib*(iunkn - 1) + [1:nlib];
                pglib(index_glib,1) = (pnods(1,:)-1)*nunkn + iunkn;
                qglib(index_glib,1) = (pnods(2,:)-1)*nunkn + iunkn;
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
    end
end
