classdef SpecimenModel < handle

    properties (Access = protected)
        geometry
    end

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.contact
                case 'thirdMedium'
                    obj = ThirdMediumSpecimenModel(cParams);
                otherwise
                    obj = FittedSpecimenModel(cParams);
            end
        end

    end

    methods (Abstract, Access = public)
        mesh = createMesh(obj,refinement,workDir,name)
        order = interpolationOrder(obj)
        p = createMaterialProperties(obj,matProp,mesh)
        s = addContact(obj,s)
        f = createEnergyEvaluator(obj,p,mesh)
        [funs,names] = functionsToPrint(obj,u)
        itDoes = checksOverlaps(obj)
    end

end
