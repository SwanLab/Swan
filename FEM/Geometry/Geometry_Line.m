classdef Geometry_Line < Geometry
    
    properties (GetAccess = public, SetAccess = private)
        dvolu
    end
    
    properties (Access = private)
        mesh
        interpolationVariable  
        quadrature      
    end
    
    methods (Access = public)
        
        function obj = Geometry_Line(cParams)
            obj.init(cParams)
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.initGeometry(interpV,quad);
            obj.computeDvolu();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
        function initGeometry(obj,interpV,quad)
            obj.interpolationVariable = interpV;
            obj.quadrature = quad;
            obj.computeShapeFunctions();
        end
        
        function computeShapeFunctions(obj)
            xpg = obj.quadrature.posgp;
            obj.interpolationVariable.computeShapeDeriv(xpg)
            obj.mesh.interpolation.computeShapeDeriv(xpg);
        end
        
        function computeDvolu(obj)
            nGaus  = obj.quadrature.ngaus;
            nDime  = obj.mesh.ndim;            
            drDtxi = zeros(nGaus,obj.mesh.nelem);
            xp     = permute(obj.mesh.coordElem,[1 3 2]);
            deriv  = obj.mesh.interpolation.deriv(1,:,:);
            dShapes = permute(deriv,[3 2 1]);
            for idime = 1:nDime
                x      = xp(:,:,idime);
                dxDtxi = dShapes*x;
                drDtxi = drDtxi + (dxDtxi).^2;
            end
            w(:,1) = obj.quadrature.weigp;
            dv =  bsxfun(@times,w,sqrt(drDtxi));
            obj.dvolu = dv';
        end        
        
    end
    
    
    
end