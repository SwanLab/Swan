classdef P0Function < FeFunction
    
    properties (Access = public)
       fElem
    end

    properties (Access = private)
        % ...
    end
    
    methods (Access = public)
        
        function obj = P0Function(cParams)
            obj.init(cParams);
            obj.createFvaluesByElem();
        end

        function fxV = interpolateFunction(obj, xV)
            % Its a p0 function, so no true need to interpolate -- the
            % value is constant
        end

        function plot(obj, m)
            [mD, fD] = obj.createDiscontinuousP0(m);
            coord  = mD.coord;
            connec = mD.connec;
            figure()
            trisurf(connec, coord(:,1), coord(:,2), fD)
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fElem = cParams.fElem;
        end

        function [mD, fD] = createDiscontinuousP0(obj, m)
            dim = 1;
            fEl = squeeze(obj.fElem(dim,:,:));
            mD = m.createDiscontinousMesh();
            nnodeElem = mD.nnodeElem;
            fRepeated = zeros(size(fEl,1), nnodeElem);
            for iNode = 1:nnodeElem
                fRepeated(:,iNode) = fEl;
            end
            fD = transpose(fRepeated);
        end

        function createFvaluesByElem(obj)
            f = obj.fElem;
            nElem = size(f,1);
            nDime = size(f,2);
            obj.fElem = reshape(f',[nDime, 1, nElem]);
        end

    end

end