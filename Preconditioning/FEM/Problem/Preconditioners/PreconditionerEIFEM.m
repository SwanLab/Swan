classdef PreconditionerEIFEM < handle

    properties (GetAccess = public, SetAccess = private)
        Fext
        EIFEMsolver
    end

    properties (Access = private)
        EIFEMfilename
        weight
        
        
    end

    properties (Access = private)
        LHS
        bcApplier
        dir
        ddDofManager
        dMesh
        iter
    end

    methods (Access = public)

        function obj = PreconditionerEIFEM(cParams)
            obj.init(cParams);
        end

        function [z,uCoarse] = apply(obj,r,uk)
% %             uk = obj.bcApplier.reducedToFullVectorDirichlet(uk);
% %             uk = obj.ddDofManager.global2local(uk);  %dissemble
% %             uk = reshape(uk,[],1);
%             Rd = obj.computeDiscontinousField(r);
%             uD = obj.EIFEMsolver.apply(Rd);
% %             u = reshape(uD,[],1);
% %             EIFEMtesting.plotSolution(u+uk,obj.dMesh,21,5,obj.iter,[],0)
%             obj.iter = obj.iter+1;
%             uC = obj.computeContinousField(uD);
%             z  = uC; 
              obj.Fext = r;
              [z,uCoarse] = obj.solve();
        end

        function [uC,uCoarse]= solve(obj)
            Rd           = obj.computeDiscontinousField(obj.Fext);
            [uD,uCoarse] = obj.EIFEMsolver.apply(Rd);
            uC           = obj.computeContinousField(uD);
        end


        function updateDownscaling(obj,mu)
            obj.EIFEMsolver.updateDownscaling(mu)
        end

        function computeLHS(obj,mu)
           obj.EIFEMsolver.computeLHS(mu)
        end
        
        function dK = computeGradK(obj,mu)
           dK = obj.EIFEMsolver.computeGradK(mu);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ddDofManager = cParams.ddDofManager;
            obj.EIFEMsolver  = cParams.EIFEMsolver;
            obj.bcApplier    = cParams.bcApplier;
            obj.weight       = 0.5;
            obj.dMesh        = cParams.dMesh;
            obj.iter         = 1;
            if isfield(cParams,'Fext')
                obj.Fext         = cParams.Fext;
            end
        end


        function Rd = computeDiscontinousField(obj,Rc)
            RG = obj.bcApplier.reducedToFullVectorDirichlet(Rc);
            Rd = obj.ddDofManager.global2local(RG);  %dissemble
            Rd = obj.ddDofManager.scaleInterfaceValues(Rd,obj.weight);    %scale
        end

        function uC = computeContinousField(obj,uD)
            fS  = obj.ddDofManager.scaleInterfaceValues(uD,obj.weight);         %scale
            uC  = obj.ddDofManager.AssembleLocal2GlobalVector(fS);
%             fG  = obj.ddDofManager.local2global(fS);   %assemble
%             uC  = sum(fG,2);                           %assemble
            uC  = obj.bcApplier.fullToReducedVectorDirichlet(uC);
        end

    end

end