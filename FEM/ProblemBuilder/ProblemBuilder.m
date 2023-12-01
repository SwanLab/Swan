classdef ProblemBuilder < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        type
        stiffness
        forces
        boundaryConditions
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ProblemBuilder(cParams)
            obj.init(cParams)
        end

        function [LHS, RHS] = compute(obj)
            lhsC = obj.createLHSconditions();
            LHS  = obj.assembleLHS(lhsC);
            rhsC = obj.createRHSconditions();
            RHS  = obj.assembleRHS(rhsC);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type               = 'MONOLITHIC';
            obj.stiffness          = cParams.stiffness;
            obj.forces             = cParams.forces;
            obj.boundaryConditions = cParams.boundaryConditions;
        end

        function Ct = createLHSconditions(obj)
            dirich = obj.boundaryConditions.dirichletFun;
            % dir_dofs = sort(dirich.getDofs());
            dir_dofs = dirich.getDofs();
            nDofs = dirich.fun.nDofs;
            nDirich = length(dir_dofs);
            Ct = full(sparse(1:nDirich, dir_dofs, 1, nDirich, nDofs));
        end

        function LHS = assembleLHS(obj, Ct)
            C   = Ct';
            nC  = size(Ct,1);
            Z   = zeros(nC);
            Km  = obj.stiffness;
            LHS = [Km C; C' Z];
        end

        function Ct = createRHSconditions(obj)
            dirich = obj.boundaryConditions.dirichletFun;
            % dir_dofs = dirich.getDofs();
            dir_vals = dirich.getValues();
            % mat = [dir_dofs, dir_vals];
            % Ct        = mat(:,2); % uD
            Ct = dir_vals;
        end

        function RHS = assembleRHS(obj, Ct)
            RHS = [obj.forces; Ct];
        end

        
    end
    
end