classdef ShFunWithElasticPdes < ShapeFunctional
    
    properties (Access = protected)
        interpolation
        physicalProblem
    end
    
    properties (Access = private)
        orientationUpdater
        alpha
    end
    
    properties (Access = public)
        regDesignVariable
    end
    
    methods (Access = public)
        
        function computeFunctionAndGradient(obj)
            obj.updateAlpha();
            obj.computeFunction();
            obj.updateAlpha();
            obj.computeGradient();
            obj.updateAlpha();
        end
        
        function computeFunction(obj)
            obj.updateHomogenizedMaterialProperties();
            obj.solveState();
            obj.computeFunctionValue();
            obj.normalizeFunction();
        end
        
        function computeGradient(obj)
            obj.solveAdjoint();
            obj.computeGradientValue();
            obj.filterGradient();
            obj.normalizeGradient();
        end
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
        end
        
        function f = getRegularizedDesignVariable(obj)
            f = obj.regDesignVariable{1:end-1};
        end
        
        function q = getQuad(obj)
            q = obj.physicalProblem.getQuadrature();
        end
        
    end
    
    methods (Access = protected)
        
        function createOrientationUpdater(obj)
            cParams.type = 'MinimumEigenValue';
            obj.orientationUpdater = OrientationUpdater.create(cParams);
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
        end
        
        function filterDesignVariable(obj)
            nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
            x  = obj.designVariable.value;
            xf = cell(obj.designVariable.nVariables,1);
            for ivar = 1:obj.designVariable.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xf{ivar} = obj.filter.getP0fromP1(xs);
            end
            xf{ivar+1} = obj.designVariable.alpha;
            obj.regDesignVariable = xf;
        end
        
        function filterGradient(obj)
            g = obj.gradient;
            gf = zeros(size(obj.Msmooth,1),obj.nVariables);
            for ivar = 1:obj.nVariables
                gs = g(:,:,ivar);
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            gf = obj.Msmooth*gf;
            g = gf(:);
            obj.gradient = g;
        end
        
        function v = getPdeVariableToPrint(obj,p)
            cParams.physicalProblem = p;
            g = PdeVariableToPrintGetter(cParams);
            v = g.compute();
        end
        
        function fP = addHomogPrintVariablesNames(obj,fP)
            fH = obj.homogenizedVariablesComputer.createPrintVariables();
            nP = numel(fP);
            for i = 1:numel(fH)
                fP{nP+i} = fH{i};
            end
        end
      
      function updateAlpha(obj)
          if isequal(obj.designVariable.type,'MicroParams')
            if isfield(obj.physicalProblem.variables,'principalStress')
            cParams.pD = obj.physicalProblem.variables.principalDirections;
            cParams.pS = obj.physicalProblem.variables.principalStress;
            obj.orientationUpdater.compute(cParams);
            alpha = obj.orientationUpdater.alpha;
            obj.designVariable.alpha = alpha;
            end
          end
      end
        
      function fP = addHomogVariables(obj,fP)
          fH = obj.homogenizedVariablesComputer.addPrintableVariables(obj.designVariable);
          for i = 1:numel(fH)
              fP{end+1} = fH{i};
          end
      end
        
    end
    
    methods (Access = private)

        function s = createFEMparameters(obj, file)
            gidParams = obj.createGiDparameters(file);
            s.dim       = gidParams.pdim;
            s.type      = gidParams.ptype;
            s.scale     = gidParams.scale;
            s.mesh      = gidParams.mesh;
            s.dirichlet = gidParams.dirichlet;
            s.pointload = gidParams.pointload;
        end

        function gidParams = createGiDparameters(obj,file)
            gidReader = FemInputReader_GiD();
            gidParams = gidReader.read(file);
        end
    
    end

    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        computeGradientValue(obj)
        solveState(obj)
        solveAdjoint(obj)
    end
end