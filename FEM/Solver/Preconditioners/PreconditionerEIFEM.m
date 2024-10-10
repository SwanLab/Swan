classdef PreconditionerEIFEM < handle

    properties (Access = public)

    end

    properties (Access = private)
        EIFEMfilename
        weight
        EIFEMsolver
    end

    properties (Access = private)
        coarseMesh
        LHS
        bcApplier
        dir
        ddDofManager
        nSubdomains
    end

    methods (Access = public)

        function obj = PreconditionerEIFEM(cParams)
            obj.init(cParams);
            obj.createEIFEM();
        end

        function z = solveEIFEM(obj,r)
            RGsbd = obj.computeSubdomainResidual(r);
            uSbd =  obj.EIFEMsolver.apply(RGsbd);
            z = obj.computeContinousField(uSbd);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.LHS          = cParams.LHS;
            obj.ddDofManager = cParams.ddDofManager;
            obj.coarseMesh   = cParams.coarseMesh;
            obj.bcApplier     = cParams.bcApplier;
            obj.dir = cParams.dir;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            obj.EIFEMfilename = 'DEF_Q4porL_1.mat';
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            obj.weight       = 0.5;
        end

        function createEIFEM(obj)
            filename        = obj.EIFEMfilename;
            RVE             = TrainedRVE(filename);
            s.RVE           = RVE;
            s.mesh          = obj.coarseMesh;
            s.DirCond       = obj.dir;
            eifem           = EIFEM(s);
            obj.EIFEMsolver = eifem;
        end

        function Rd = computeSubdomainResidual(obj,Rc)
            RG    = obj.bcApplier.reducedToFullVectorDirichlet(Rc);
            Rd = obj.ddDofManager.global2local(RG);  %dissemble
            Rd = obj.ddDofManager.scaleInterfaceValues(Rd,obj.weight);    %scale
        end

        function uC = computeContinousField(obj,uD)
            fS  = obj.ddDofManager.scaleInterfaceValues(uD,obj.weight);         %scale
            fG  = obj.ddDofManager.local2global(fS);   %assemble
            uC  = sum(fG,2);                           %assemble
            uC  = obj.bcApplier.fullToReducedVectorDirichlet(uC);
        end

    end

end