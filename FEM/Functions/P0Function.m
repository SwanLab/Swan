classdef P0Function < FeFunction
    
    properties (Access = public)
    end

    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = P0Function(cParams)
            obj.init(cParams);
            obj.createFvaluesByElem();
        end

        function fxV = evaluate(obj, xV)
            % Its a p0 function, so no true need to interpolate -- the
            % value is constant
            fxV = obj.fValues;
        end

        function fD = computeP1DiscontinuousFunction(obj)
            dim = 1;
            fEl = squeeze(obj.fValues(dim,:,:));
            mD = m.createDiscontinousMesh();
            nnodeElem = mD.nnodeElem;
            fRepeated = zeros(size(fEl,1), nnodeElem);
            for iNode = 1:nnodeElem
                fRepeated(:,iNode) = fEl;
            end
            fD = transpose(fRepeated);
        end

        function plot(obj, m)
%             f1 = obj.computeP1DiscontinuousFunction(); %public
%             f1.plot();
            [mD, fD] = obj.createDiscontinuousP0(m);
            coord  = mD.coord;
            connec = mD.connec;
            figure()
            trisurf(connec, coord(:,1), coord(:,2), fD)
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fValues = cParams.fValues;
        end

        function [mD, fD] = createDiscontinuousP0(obj, m)
            dim = 1;
            fEl = squeeze(obj.fValues(dim,:,:));
            mD = m.createDiscontinuousMesh();
            nnodeElem = mD.nnodeElem;
            fRepeated = zeros(size(fEl,1), nnodeElem);
            for iNode = 1:nnodeElem
                fRepeated(:,iNode) = fEl;
            end
            fD = transpose(fRepeated);
        end

        function createFvaluesByElem(obj)
            f = obj.fValues;
            nElem = size(f,1);
            nDime = size(f,2);
            obj.fValues = reshape(f',[nDime, 1, nElem]);
        end

    end

end