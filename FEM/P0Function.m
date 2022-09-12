classdef P0Function < FeFunction
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        quadrature
        dim
    end
    
    properties (Access = private)
       mesh 
       fValues
       quadOrder
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
        
        function plot(obj, m, f)
%             disMesh = m.createDiscontinousMesh();
%             p0funct = obj.preparePlottingFunction(disMesh, f);
%             coor = disMesh.coord;
%             conn = disMesh.connecM;
%             trisurf(conn, coor(:,1), coor(:,2), p0funct)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.fValues   = cParams.fValues;
            obj.quadOrder = 'LINEAR';
%             obj.connec    = cParams.connec;
        end

        function p0F = preparePlottingFunction(obj, dM, f)
            nnodeElem = dM.nnodeElem;
            fRepeated = zeros(size(f,1), nnodeElem);
            for iNode = 1:nnodeElem
                fRepeated(:,iNode) = f;
            end
            fRepeated = transpose(fRepeated);
            p0F = fRepeated(:);
        end

        function createFvaluesByElem(obj)
        end

    end
    
end