classdef ThirdMediumSpecimenModel < SpecimenModel

    properties (Access = private)
        regularization
        mediumStiffness
        density
    end

    methods (Access = public)

        function obj = ThirdMediumSpecimenModel(cParams)
            obj.geometry        = cParams.geometry;
            obj.regularization  = cParams.regularization;
            obj.mediumStiffness = cParams.mediumStiffness;
        end

        function mesh = createMesh(obj,refinement,workDir,name)
            mesh = obj.geometry.createBackgroundMesh(refinement,workDir,name);
            obj.density = obj.geometry.computeDensity(mesh);
        end

        function order = interpolationOrder(~)
            order = 'P2';
        end

        function p = createMaterialProperties(obj,matProp,~)
            p.type           = 'ThirdMedium';
            p.mu             = matProp.mu;
            p.lambda         = matProp.lambda;
            p.density        = obj.density;
            p.initialGamma   = 1e-2;
            p.finalGamma     = obj.mediumStiffness;
            p.gammaReduction = 10;
            p.regularization = obj.regularization;
        end

        function s = addContact(~,s)
        end

        function f = createEnergyEvaluator(~,p,mesh)
            s.mesh     = mesh;
            s.material = p;
            s.material.initialGamma = p.finalGamma;
            f = ThirdMediumFunctional(s);
        end

        function [funs,names] = functionsToPrint(obj,u)
            funs  = {u,obj.density};
            names = {'fValues','density'};
        end

        function itDoes = checksOverlaps(~)
            itDoes = false;
        end

    end

end
