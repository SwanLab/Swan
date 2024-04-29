classdef Multigrid < handle


    properties (Access = private)
        nDimf
        data
        nLevel
        mesh
        material
        boundaryConditions
        BcApplierFine 
        tol
        type
        scale  
        pdim
        RHS
        LHS
        coarseMeshes
        interpolator
%         coarseBc
        bcApplierCoarse
        coarseLHS
        coarseRHS
        solver
        solverCase
    end

    methods (Access = public)

        function obj = Multigrid(cParams)
            obj.init(cParams);            
            obj.createFEMlevel();
        end

        function u = solve(obj)
            obj.bcApplierCoarse(obj.nLevel+1) = obj.BcApplierFine;
            obj.coarseLHS{1,obj.nLevel + 1}   = obj.LHS;
            obj.coarseRHS{1,obj.nLevel + 1}   = obj.RHS;
            u                                 = 0*obj.coarseRHS{obj.nLevel + 1};
            level                             = obj.nLevel + 1;
            b                                 = obj.RHS;
            Res                               = 1;
            
            while norm(Res, inf) >= obj.tol
                u = obj.vCycle(u,b,level);
                Res = obj.RHS - obj.LHS * u;
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
            obj.BcApplierFine      = cParams.bcApplier;
            obj.material           = cParams.material;
            obj.LHS                = cParams.LHS;
            obj.RHS                = cParams.RHS;
            obj.type               = cParams.type ;%'ELASTIC';
            obj.scale              = cParams.scale; %'MACRO';
            obj.pdim               = cParams.dim; %'2D';
            obj.solverCase         = cParams.solverCase;
            obj.nDimf              = cParams.nDimf;
        end

        function createFEMlevel(obj)
            s.coarseMeshes       = obj.coarseMeshes;
            s.nDimf              = obj.nDimf;
            s.nLevel             = obj.nLevel;
            s.type               = obj.type;
            s.scale              = obj.scale;
            s.pdim               = obj.pdim;
            s.tol                = obj.tol;
            s.solverCase         = obj.solverCase;
            FEM                  = FemCreator(s);
            obj.bcApplierCoarse  = FEM.bcApplier;
            obj.coarseLHS        = FEM.LHS;
            obj.coarseRHS        = FEM.RHS;
            obj.solver           = FEM.solver;
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
                bcApplierold = obj.bcApplierCoarse(level-1);
                bcApplier = obj.bcApplierCoarse(level);
                u      = obj.solver{level}.solve(LHS,RHS,u);
                r      = b - LHS*u;
                Rr     = obj.interpolate(r,bcApplierold,bcApplier,intOld);
                ur     = obj.interpolate(u,bcApplierold,bcApplier,intOld);
                er     = obj.vCycle(0*ur, Rr, level - 1);
                e      = obj.restriction(er,bcApplierold,bcApplier,intOld);
                u      = u + e; 
                u      = obj.solver{level}.solve(LHS,RHS,u);
             end

        end

        function fF = restriction(obj,fC,bcApplierold,bcApplier,IOld)
                fF = bcApplierold.reducedToFullVectorDirichlet(fC);
                fF = reshape(fF,obj.nDimf,[])';
                fF = IOld * fF; 
                fF = reshape(fF',[],1);
                fF = bcApplier.fullToReducedVectorDirichlet(fF);
        end

        function fCoarse = interpolate(obj,fFine,bcApplierold,bcApplier,IOld)
                fFine   = bcApplier.reducedToFullVectorDirichlet(fFine);
                fFine   = reshape(fFine,obj.nDimf,[])';
                fCoarse = IOld' * fFine;
                fCoarse = reshape(fCoarse',[],1);
                fCoarse = bcApplierold.fullToReducedVectorDirichlet(fCoarse);
        end
        
         function plotRes(obj,res,mesh,bcApplier,numItr)
            xFull = bcApplier.reducedToFullVectorDirichlet(res);
            s.fValues = reshape(xFull,obj.nDimf,[])';
            s.mesh = mesh;
            if obj.nDimf<3
                s.fValues(:,end+1) = 0;
            end
             s.order   = 'P1';
%             s.ndimf = obj.nDimf;
            uFeFun = LagrangianFunction(s);
%             xF = P1Function(s);
            %xF.plot();
            uFeFun.print('uPrueva','Paraview')
            fclose('all');
        end
        
    end
end

