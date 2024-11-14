classdef MultiMaterialCharacteristicFunction < handle

    properties (Access = private)
        designVariable
        mesh
    end

    methods (Access = public)
        function obj = MultiMaterialCharacteristicFunction(cParams)
            obj.init(cParams);
        end

        function charFun = compute(obj)
            nLS     = length(obj.designVariable.levelSets);
            chiLoc  = obj.computeParticularCharFunValues();
            chi     = obj.computeCharFunValues(chiLoc);
            charFun = cell(nLS+1,1);
            for i = 1:nLS+1
                charFun{i} = LagrangianFunction.create(obj.mesh,1,'P0');
                charFun{i}.fValues = chi(:,i);
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.designVariable = cParams.designVariable;
        end

        function chiProj = computeParticularCharFunValues(obj)
            nLS     = length(obj.designVariable.levelSets);
            nElem   = obj.mesh.nelem;
            chiProj = zeros(nElem,nLS);
            for i = 1:nLS
                x            = obj.designVariable.levelSets{i};
                chi          = x.obtainDomainFunction();
                chiV         = chi.project('P0');
                chiProj(:,i) = chiV.fValues;
            end
        end

        function phi = computeCharFunValues(obj,chi)
            nLS = length(obj.designVariable.levelSets);
            phi = zeros(size(chi,1),nLS+1);
            for i=1:nLS-1
                phi(:,i) = (1 - chi(:,i+1)).*prod(chi(:,1:i),2);
            end
            phi(:,nLS)   = prod(chi,2);
            phi(:,nLS+1) = (1 - chi(:,1));
        end
    end
end