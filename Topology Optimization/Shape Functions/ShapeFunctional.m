classdef ShapeFunctional < handle
    
    properties (Access = public)
        value
        gradient
        filter
        Msmooth
        dvolu
        value0
    end
    
    properties (Access = protected)
        homogenizedVariablesComputer
        nVariables
        designVariable
        target_parameters;
    end

    properties (Access = private)
        dim
        mesh
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = ShapeFunctional_Factory();
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
            obj.createFilter(cParams);
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
        
        function createFilter(obj,cParams)
            s = cParams.filterParams;
            s.femSettings.mesh = s.mesh;
            s.designVariable = cParams.designVariable;
            obj.filter = NewFilter.create(s);
            obj.filter.preProcess();
        end
        
        function createMsmoothAndDvolu(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.computeDimensions();
            q = Quadrature.set(cParams.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.Msmooth = obj.computeMassMatrix();
            obj.dvolu = cParams.mesh.computeDvolume(q)';
        end

        function computeDimensions(obj)
            s.type = 'Scalar';
            s.name = 'x';
            s.mesh = obj.mesh;
            dims   = DimensionVariables.create(s);
            obj.dim = dims;
        end
        
        function M = computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            M = LHS.compute();
        end

    end
    
end