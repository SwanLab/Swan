classdef testUnfittedIntegration_ExternalIntegrator_Composite < testUnfittedIntegration_ExternalIntegrator
    
    methods (Access = protected)
        
        function B = sumResults(obj,A)
            Bi = obj.sumInteriorResults(A);
            Bb = obj.sumBoxFacesResults(A);
            B = Bi + Bb;
        end
        
    end
    
    methods (Access = private)
        
        function B = sumBoxFacesResults(obj,A)
            B = 0;
            for iactive = 1:obj.mesh.nActiveBoxFaces
                iface = obj.mesh.activeBoxFaceMeshesList(iactive);
                B = B + sum(A.boxFacesIntegrals{iface});
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function B = sumInteriorResults(A)
            B = sum(A.interiorIntegral);
        end
        
    end
    
end

