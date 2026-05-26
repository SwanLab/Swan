classdef VolumeConstraintMicroDensity < handle

    properties (Access = private)
        mesh
        volumeTarget
        volume
    end
    
    methods (Access = public)
        function obj = VolumeConstraintMicroDensity(cParams)
            obj.init(cParams);
        end
        
        function [J, dJ] = computeFunctionAndGradient(obj, x)
            [V, dV] = obj.volume.computeFunctionAndGradient(x);
            J = obj.computeFunction(V);
            
            % dV{1} é o gradiente em relação a ρ (supondo que VolumeFunctionalMicroDensity já retorna correctamente)
            dJ_rho = obj.computeGradient(dV{2});
            
            % Gradiente em b = ZERO (volume não depende de b)
            dJ_b = dJ_rho.copy();
            dJ_b.setFValues(zeros(size(dJ_rho.fValues)));
            
            % Concatenar na ordem correcta: {b, ρ}
            dJ{1} = dJ_b;   % gradiente em b
            dJ{2} = dJ_rho; % gradiente em ρ
        end  
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.volume       = VolumeFunctionalMicroDensity(cParams);
        end

        function J = computeFunction(obj, V)
            vTar = obj.volumeTarget;
            J = V / vTar - 1;
        end

        function dJ = computeGradient(obj, dV)
            vTar = obj.volumeTarget;
            dJ   = dV.copy();
            dJ.setFValues(dV.fValues / vTar);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end