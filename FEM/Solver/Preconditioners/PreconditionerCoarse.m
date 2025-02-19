classdef PreconditionerCoarse < handle

    properties (Access = public)

    end

    properties (Access = private)
        EIFEMfilename
        weight
        Coarsesolver
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

        function obj = PreconditionerCoarse(cParams)
            obj.init(cParams);
        end

        function z = apply(obj,r)
            Rd = obj.computeDiscontinousField(r);
            uD = obj.Coarsesolver.apply(Rd);
            uC = obj.computeContinousField(uD);
            z  = uC; 
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ddDofManager = cParams.ddDofManager;
            obj.Coarsesolver  = cParams.Coarsesolver;
            obj.bcApplier    = cParams.bcApplier;
            obj.weight       = 0.5;
            obj.dMesh        = cParams.dMesh;
            obj.iter         = 1;
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