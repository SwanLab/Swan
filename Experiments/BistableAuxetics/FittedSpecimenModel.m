classdef FittedSpecimenModel < SpecimenModel

    properties (Access = private)
        contact
    end

    methods (Access = public)

        function obj = FittedSpecimenModel(cParams)
            obj.geometry = cParams.geometry;
            obj.contact  = cParams.contact;
        end

        function mesh = createMesh(obj,refinement,workDir,name)
            mesh = obj.geometry.createMesh(refinement,workDir,name);
        end

        function order = interpolationOrder(~)
            order = 'P1';
        end

        function p = createMaterialProperties(~,matProp,mesh)
            p.type   = 'Neohookean';
            p.mu     = ConstantFunction.create(matProp.mu,mesh);
            p.lambda = ConstantFunction.create(matProp.lambda,mesh);
        end

        function s = addContact(obj,s)
            w = obj.geometry.hingeWidth;
            switch obj.contact
                case 'penalty'
                    s.contact.penalty = 1e3;
                    s.contact.depth   = w;
                case 'augmentedLagrangian'
                    s.contact.penalty = 1e2;
                    s.contact.depth   = w;
                    s.contact.penetrationTolerance = 0.01*w;
                case 'barrier'
                    s.contact           = obj.geometry.computeBarrierParameters();
                    s.contact.stiffness = 1;
                otherwise
                    return
            end
            s.contact.type = obj.contact;
        end

        function f = createEnergyEvaluator(~,p,mesh)
            s.mesh     = mesh;
            s.material = p;
            f = NeohookeanFunctional(s);
        end

        function [funs,names] = functionsToPrint(~,u)
            funs  = {u};
            names = {'fValues'};
        end

        function itDoes = checksOverlaps(~)
            itDoes = true;
        end

    end

end
