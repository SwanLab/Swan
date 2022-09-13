classdef P0Function < FeFunction
    
    properties (Access = private)
        meshDisc
        fDisc
    end
    
    properties (Access = private)
       mesh 
       fNodes
    end
    
    methods (Access = public)
        
        function obj = P0Function(cParams)
            obj.init(cParams);
            obj.createDiscontinuousP0();
        end

        function plot(obj)
            coord  = obj.meshDisc.coord;
            connec = obj.meshDisc.connec;
            figure()
            trisurf(connec, coord(:,1), coord(:,2), obj.fDisc)
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.fNodes = cParams.fNodes;
        end

        function createDiscontinuousP0(obj)
            obj.meshDisc = obj.mesh.createDiscontinousMesh();
            nnodeElem = obj.meshDisc.nnodeElem;
            fRepeated = zeros(size(obj.fNodes,1), nnodeElem);
            for iNode = 1:nnodeElem
                fRepeated(:,iNode) = obj.fNodes;
            end
            obj.fDisc = transpose(fRepeated);
        end

    end

end