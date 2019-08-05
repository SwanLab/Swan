classdef testUnfittedIntegration_ExternalIntegrator_Composite < testUnfittedIntegration_ExternalIntegrator
    
    methods (Access = protected)
        
        function intF = sumResults(obj,f)
%                         Bi = obj.sumInteriorResults(f);
%                         Bb = obj.sumBoxFacesResults(f);
            %             B = Bi + Bb;
            intF = 0;
            for iInt = 1:obj.integrator.nInt
                intF = intF + sum(f{iInt});
            end
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

