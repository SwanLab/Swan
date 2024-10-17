classdef PreconditionerEIFEM < handle

    properties (Access = public)

    end

    properties (Access = private)
        EIFEMfilename
        weight
        EIFEMsolver
    end

    properties (Access = private)
        LHS
        bcApplier
        dir
        ddDofManager
    end

    methods (Access = public)

        function obj = PreconditionerEIFEM(cParams)
            obj.init(cParams);
        end

        function z = apply(obj,r)
            Rd = obj.computeDiscontinousField(r);
            uD = obj.EIFEMsolver.apply(Rd);
            uC = obj.computeContinousField(uD);
            z  = uC; 
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ddDofManager = cParams.ddDofManager;
            obj.EIFEMsolver  = cParams.EIFEMsolver;
            obj.bcApplier    = cParams.bcApplier;
            obj.weight       = 0.5;
        end


        function Rd = computeDiscontinousField(obj,Rc)
            RG = obj.bcApplier.reducedToFullVectorDirichlet(Rc);
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