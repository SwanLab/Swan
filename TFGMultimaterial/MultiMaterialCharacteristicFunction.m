classdef MultiMaterialCharacteristicFunction < handle

    properties (Access = private)
        designVariable
        mesh
    end

    methods (Access = public)
        function obj = MultiMaterialCharacteristicFunction(cParams)
            obj.init(cParams);
        end


        function [fi, tfi] = computeAtNodesAndElements(obj)
            fi  = obj.computeFunction('P1');
            tfi = obj.computeFunction('P0');            
        end
    end
    

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.designVariable = cParams.designVariable;
        end

        function charFun = computeFunction(obj,order)
            nLS     = length(obj.designVariable.levelSets);
            chiLoc  = obj.computeParticularCharFunValues(order);
            chi     = obj.computeCharFunValues(chiLoc);
            charFun = LagrangianFunction.create(obj.mesh,nLS+1,'P1');
            charFun.fValues = chi;
        end

        function chiProjectedVals = computeParticularCharFunValues(obj,order)
            nLS              = length(obj.designVariable.levelSets);
            chiProjectedVals = [];
            for i = 1:nLS
                x   = obj.designVariable.levelSets{i};
                chi = x.obtainDomainFunction();
                chi = chi.project(order);
                chiProjectedVals = [chiProjectedVals,chi.fValues];
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