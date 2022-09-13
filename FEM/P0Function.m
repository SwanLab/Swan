classdef P0Function < FeFunction
    
    properties (Access = private)
        meshDisc
        fDisc
    end
    
    properties (Access = private)
       mesh 
       fElem
       fByElem
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

        function plot(obj)
            obj.createDiscontinuousP0();
            coord  = obj.meshDisc.coord;
            connec = obj.meshDisc.connec;
            figure()
            trisurf(connec, coord(:,1), coord(:,2), obj.fDisc)
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.fElem = cParams.fElem;
        end

        function createDiscontinuousP0(obj)
            dim = 1;
            fEl = squeeze(obj.fElem(dim,:,:));
            obj.meshDisc = obj.mesh.createDiscontinousMesh();
            nnodeElem = obj.meshDisc.nnodeElem;
            fRepeated = zeros(size(fEl,1), nnodeElem);
            for iNode = 1:nnodeElem
                fRepeated(:,iNode) = fEl;
            end
            obj.fDisc = transpose(fRepeated);
        end

        function createFvaluesByElem(obj)
            f = obj.fElem;
            nElem = size(f,1);
            nDime = size(f,2);
            obj.fElem = reshape(f',[nDime, 1, nElem]);
        end

    end

end