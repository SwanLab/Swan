classdef Multigrid < handle


    properties (Access = private)
        nDimf
        data
        nLevel
        mesh
        material
        boundaryConditions
        tol
        type
        scale  
        pdim
        RHS
        LHS
        coarseMeshes
        interpolator
        coarseBc
        coarseLHS
        coarseRHS
        solver
    end

    methods (Access = public)

        function obj = Multigrid(cParams)
            obj.init(cParams);            
            obj.createFEMlevel();
        end

        function u = solve(obj)
            obj.coarseBc{1,obj.nLevel + 1}  = obj.boundaryConditions;
            obj.coarseLHS{1,obj.nLevel + 1} = obj.LHS;
            obj.coarseRHS{1,obj.nLevel + 1} = obj.RHS;
            u                               = 0*obj.coarseRHS{obj.nLevel + 1};
            level                           = obj.nLevel + 1;
            b                               = obj.RHS;
            
            while norm(obj.coarseRHS{obj.nLevel + 1} - obj.coarseLHS{obj.nLevel + 1} * u, inf) >= obj.tol
                u = obj.vCycle(u,b,level);
            end

        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.coarseMeshes       = cParams.coarseMeshes;
            obj.interpolator       = cParams.interpolator;
            obj.tol                = cParams.tol;
            obj.nLevel             = cParams.nLevel;
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.bc;
            obj.material           = cParams.material;
            obj.LHS                = cParams.LHS;
            obj.RHS                = cParams.RHS;
            obj.type               = cParams.type ;%'ELASTIC';
            obj.scale              = cParams.scale; %'MACRO';
            obj.pdim               = cParams.dim; %'2D';
            obj.nDimf              = cParams.nDimf;
        end

        function createFEMlevel(obj)
            s.coarseMeshes = obj.coarseMeshes;
            s.nDimf        = obj.nDimf;
            s.nLevel       = obj.nLevel;
            s.type         = obj.type;
            s.scale        = obj.scale;
            s.pdim         = obj.pdim;
            s.tol          = obj.tol;
            FEM            = FemCreator(s);
            obj.coarseBc   = FEM.bc;
            obj.coarseLHS  = FEM.LHS;
            obj.coarseRHS  = FEM.RHS;
            obj.solver     = FEM.solver;
        end
        
        function u = vCycle(obj,u,b,level)
            if level == 1
                LHS = obj.coarseLHS{level};
                RHS = b;
                u   = obj.solver{level}.solve(LHS,RHS,u);
            else
                LHS    = obj.coarseLHS{level};                
                RHS    = b;
                intOld = obj.interpolator{level-1};
                bc     = obj.coarseBc{level};                
                bcOld  = obj.coarseBc{level-1};
                u      = obj.solver{level}.solve(LHS,RHS,u);
                r      = b - LHS*u;
                Rr     = obj.interpolate(r,bcOld,bc,intOld);
                ur     = obj.interpolate(u,bcOld,bc,intOld);
                er     = obj.vCycle(0*ur, Rr, level - 1);
                e      = obj.restriction(er,bcOld,bc,intOld);
                u      = u + e; 
                u      = obj.solver{level}.solve(LHS,RHS,u);
             end

        end

        function fF = restriction(obj,fC,bcOld,bc,IOld)
                fF = bcOld.reducedToFullVector(fC);
                fF = reshape(fF,2,[])';
                fF = IOld * fF; 
                fF = reshape(fF',[],1);
                fF = bc.fullToReducedVector(fF);
        end

        function fCoarse = interpolate(obj,fFine,bcOld,bc,IOld)
                fFine   = bc.reducedToFullVector(fFine);
                fFine   = reshape(fFine,2,[])';
                fCoarse = IOld' * fFine;
                fCoarse = reshape(fCoarse',[],1);
                fCoarse = bcOld.fullToReducedVector(fCoarse);
        end
        
    end
end

