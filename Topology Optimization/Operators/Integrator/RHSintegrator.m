classdef RHSintegrator < handle

    properties (Access = protected)
        mesh
        quadrature
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = RHSintegratorFactory();
            obj = f.create(s);
        end

    end
    
    methods (Access = public)
        
%         function obj = RHSintegrator(cParams)
%         end

        % function createQuadrature(obj)
        %     q = Quadrature.set(obj.mesh.type);
        %     q.computeQuadrature('LINEAR');
        %     obj.quadrature = q;
        % end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
        end
        
    end
    
end