classdef Geometry_Line < Geometry
    
    methods (Access = public)
        
        function obj = Geometry_Line(cParams)
            obj.permutation = [2 3 1];            
            obj.init(cParams);
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.initGeometry(interpV,quad);
            obj.computeDvolu();
        end
        
    end
    
    methods (Access = private)
               
        function computeDvolu(obj)
            nGaus  = obj.quadrature.ngaus;
            nDime  = obj.mesh.ndim;            
            drDtxi = zeros(nGaus,obj.mesh.nelem);
            xp = obj.coordElem;
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