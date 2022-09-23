classdef P0Function < FeFunction
    
    properties (Access = public)
    end

    properties (Access = private)
        connec
        type
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

        function fDfun = computeP1DiscontinuousFunction(obj)
            fD = obj.createDiscontinuousFunction();
            s.fValues = fD;
            s.connec  = obj.connec;
            s.type    = obj.type;
            fDfun = P1DiscontinuousFunction(s);            
        end

        function plot(obj, m)
            p1DiscFun = obj.computeP1DiscontinuousFunction();
            p1DiscFun.plot(m);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.connec  = cParams.connec; % Needed for discontinuous
            obj.type    = cParams.type;
        end

        function createFvaluesByElem(obj)
            f = obj.fValues;
            nElem = size(f,1);
            nDime = size(f,2);
            obj.fValues = reshape(f',[nDime, 1, nElem]);
        end

        function fD = createDiscontinuousFunction(obj)
            dim = 1;
            ndim  = size(obj.fValues, 1);
            nnodeElem = size(obj.connec,2);
            fEl = squeeze(obj.fValues(dim,:,:));
            fRepeated = zeros(ndim, size(fEl,1), nnodeElem);
            for idim = 1:ndim
                fEl = squeeze(obj.fValues(idim,:,:));
                for iNode = 1:nnodeElem
                    fRepeated(idim, :,iNode) = fEl;
                end
            end
            fD = permute(fRepeated, [1 3 2]);
        end

    end

end