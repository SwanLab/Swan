classdef ComplianceFunctionalMicroDensity < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filterMicro
        filterDensity
        compliance
        material        % APENAS UM MATERIAL (combina b e ρ)
    end

    methods (Access = public)
        function obj = ComplianceFunctionalMicroDensity(cParams)
            obj.init(cParams);
            % NÃO criar filterDensity aqui - já vem de cParams
        end

        function [J, dJ] = computeFunctionAndGradient(obj, x)
            xD = x.obtainDomainFunction();
            
            % Filtrar b e ρ separadamente
            xB   = obj.filterMicroField(xD{1});    % b filtrado
            xRho = obj.filterDensityField(xD{2});  % ρ filtrado
            
            % Unificar e passar para o material único
            xUnified = {xB, xRho};
            obj.material.setDesignVariable(xUnified);
            
            [J, dJ] = obj.computeComplianceFunctionAndGradient();
        end
    end
    
    methods (Access = private)
        function init(obj, cParams)
            obj.mesh       = cParams.mesh;
            obj.filterMicro = cParams.filter;
            obj.filterDensity = cParams.filterDensity;  % vem do TutorialFirstSimultaneous
            obj.material    = cParams.material;
            obj.compliance  = cParams.complainceFromConstitutive;
            if isfield(cParams, 'value0')
                obj.value0 = cParams.value0;
            end
        end

        function xR = filterMicroField(obj, x)
            if iscell(x)
                nDesVar = length(x);
                xR = cell(nDesVar, 1);
                for i = 1:nDesVar
                    xR{i} = obj.filterMicro.compute(x{i}, 2);
                end
            else
                xR{1} = obj.filterMicro.compute(x, 2);
            end
        end

        function xR = filterDensityField(obj, x)
            if iscell(x)
                nDesVar = length(x);
                xR = cell(nDesVar, 1);
                for i = 1:nDesVar
                    xR{i} = obj.filterDensity.compute(x{i}, 2);
                end
            else
                xR{1} = obj.filterDensity.compute(x, 2);
            end
        end

        function [J, dJ] = computeComplianceFunctionAndGradient(obj)
            % Obter tensor e derivadas do material unificado
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();  % dC{1}=b, dC{2}=ρ
            
            % Gradiente em relação a b
            [J, dJ_b] = obj.compliance.computeFunctionAndGradient(C, {dC{1}});
            
            % Gradiente em relação a ρ
            [~, dJ_rho] = obj.compliance.computeFunctionAndGradient(C, {dC{2}});
            
            % Filtrar cada gradiente
            dJ_b   = obj.filterMicroField(dJ_b);
            dJ_rho = obj.filterDensityField(dJ_rho);
            
            % Concatenar na ordem correta
            dJ = [dJ_b, dJ_rho];
            
            % Normalizar
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J / obj.value0;
            dJ = obj.computeNonDimensionalGradient(dJ);
        end

        function dx = computeNonDimensionalGradient(obj, dx)
            refX = obj.value0;
            for i = 1:length(dx)
                dx{i}.setFValues(dx{i}.fValues / refX);
            end
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end