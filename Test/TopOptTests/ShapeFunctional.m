classdef ShapeFunctional < handle
    
    properties (Access = public)
        value
        gradient
        filter
        gradientFilter
        Msmooth
        dvolu
        value0
        shNumber
    end
    
    properties (Access = protected)
        homogenizedVariablesComputer
        nVariables
        designVariable
        target_parameters;
    end

    properties (Access = private)
        mesh
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = ShapeFunctionalFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function updateTargetParameters(obj)
%            if contains(class(obj.filter),'PDE')
%                obj.filter.updateEpsilon(obj.target_parameters.epsilon);
%            end
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.storeFilters(cParams.femSettings);
            obj.createMsmoothAndDvolu(cParams);
            obj.homogenizedVariablesComputer = cParams.homogVarComputer;
            obj.designVariable = cParams.designVariable;
            obj.target_parameters = cParams.targetParameters;
            obj.nVariables = obj.designVariable.nVariables;
        end
        
        function normalizeFunction(obj)
            if isempty(obj.value0)
                obj.value0 = obj.value;
            end
            obj.value = obj.value/abs(obj.value0);
        end
        
        function normalizeGradient(obj)
            obj.gradient = obj.gradient/abs(obj.value0);
        end
        
        function fP = addHomogPrintVariablesNames(obj,fP)
            fH = obj.homogenizedVariablesComputer.createPrintVariables();
            nP = numel(fP);
            for i = 1:numel(fH)
                fP{nP+i} = fH{i};
            end
        end

        function f = obtainDomainFunction(obj)
            switch obj.designVariable.type
                case 'Density'
                    f = obj.designVariable.fun;
                case 'LevelSet'
                    f = obj.designVariable.getCharacteristicFunction();
            end
        end
    end
    
    methods (Access = protected, Static)
        
        function fP = obtainPrintVariables(types,names)
            fP = cell(numel(types),1);
            for iV = 1:numel(types)
               fP{iV}.type = types{iV};
               fP{iV}.name = names{iV};
            end
        end
        
    end
    
    methods (Access = private)

        function storeFilters(obj,cParams)
            obj.filter         = cParams.designVariableFilter;
            obj.gradientFilter = cParams.gradientFilter;
        end

        function createMsmoothAndDvolu(obj,cParams)
            obj.mesh = cParams.mesh;
            q = Quadrature.set(cParams.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.Msmooth = obj.computeMassMatrix();
            obj.dvolu = cParams.mesh.computeDvolume(q)';
        end

        function M = computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.quadratureOrder = 'QUADRATICMASS';
            LHS = LHSintegrator.create(s);
            M = LHS.compute();
        end

    end
    
end