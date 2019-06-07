classdef DirichletConditionsApplier < BoundaryConditionsApplier
    
    properties (Access = private)
        nfields
        dof
    end
    
    methods (Access = public)
        
        function obj = DirichletConditionsApplier(cParams)
            obj.nfields = cParams.nfields;
            obj.dof     = cParams.dof;
        end
        
        function Ared = fullToReducedMatrix(obj,A)
            free = obj.computeGlobalFree();
            Ared = A(free,free);
        end
        
        function b_red = fullToReducedVector(obj,b)
            free = obj.computeGlobalFree();
            b_red = b(free);
        end
        
        function b = reducedToFullVector(obj,bfree)
            [dirichlet,uD,free] = obj.compute_global_dirichlet_free_uD();
            nsteps = length(bfree(1,:));
            ndof = sum(obj.dof.ndof);
            uD = repmat(uD,1,nsteps);
            
            b = zeros(ndof,nsteps);
            b(free,:) = bfree;
            if ~isempty(dirichlet)
                b(dirichlet,:) = uD;
            end
        end
        
        function [dirichlet,uD,free] = compute_global_dirichlet_free_uD(obj)
            uD = obj.computeUd();
            dirichlet = obj.computeDirichlet();
            free = obj.computeGlobalFree();
        end        
        
        
    end
    
    
    methods (Access = private)
        
    
        
        function uD = computeUd(obj)
            global_ndof=0;
            for ifield=1:obj.nfields
                uD{ifield,1} = obj.dof.dirichlet_values{ifield};
                global_ndof=global_ndof+obj.dof.ndof(ifield);
            end
            uD = cell2mat(uD);
        end
        
        function dirichlet = computeDirichlet(obj)
            global_ndof=0;
            for ifield=1:obj.nfields
                dirichlet{ifield,1} = obj.dof.dirichlet{ifield}+global_ndof;
                global_ndof=global_ndof+obj.dof.ndof(ifield);
            end
            dirichlet = cell2mat(dirichlet);
        end
        
        function free = computeGlobalFree(obj)
            global_ndof=0;
            free = cell(obj.nfields,1);
            for ifield=1:obj.nfields
                free{ifield,1} = obj.dof.free{ifield}' + global_ndof;
                global_ndof=global_ndof+obj.dof.ndof(ifield);
            end
            free = cell2mat(free);
        end
        
        
    end
    
end

