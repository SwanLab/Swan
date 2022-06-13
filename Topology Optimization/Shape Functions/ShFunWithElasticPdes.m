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
        
        function createEquilibriumProblem(obj,fileName)
            a.fileName = fileName;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
            obj.physicalProblem = FEM.create(s);
            obj.initPrincipalDirections();
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
        end
        
        function filterDesignVariable(obj)
            switch obj.designVariable.type
                case 'DensityEigModes'
                    x  = obj.designVariable.getDensity();
                    nx = length(x)/obj.designVariable.nVariables;
                otherwise
                    nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
                    x  = obj.designVariable.value;
            end
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
            obj.printInGiD(g);
            switch obj.designVariable.type
                case 'DensityEigModes'
                    obj.gradient(end+1,1) = 1;
                otherwise
                    
            end
%             e1 = obj.filterEp(1);
%             e2 = obj.filterEp(2);
%             e3 = obj.filterEp(3);
%             e = zeros(length(e1),3);
%             e(:,1) = e1;
%             e(:,2) = e2;
%             e(:,3) = e3;
%             obj.ePrint = e;
%             obj.printInGiD();

%         g = filterEp(obj,obj.v,1);
%         obj.printInGiD(g);


        end

        function g = filterEp(obj,E,num)
            g = E(num,:);
            gf = zeros(size(obj.Msmooth,1),3);
            for ivar = 1:1
                gs = g(:,:,ivar);
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            gf = obj.Msmooth*gf;
            g = gf(:);
            
        end

        function printInGiD(obj,g)
            fileName = 'EigModesG2';
            m = obj.mesh;
            quad = Quadrature.set(m.type);
            quad.computeQuadrature('LINEAR');
            dI = obj.createPostProcessDataBase(m,fileName);
            dI.ndim   = 2;
            dI.pdim   = '2D';
            dI.ptype  = 'ELASTIC';
            dI.name = '';
            dI.printers = 'DensityGauss'; % 'HomogenizedTensor'
            p = Postprocess('VectorField',dI);
            % f.u = obj.mode1;
            zerosss = zeros(length(g),1);
            dP.fields.u = [obj.gradient zerosss]; %obj.ePrint; %[obj.gradient zerosss];
            dP.quad   = quad;
            iter = 0;
            p.print(iter,dP);
        end

        function d = createPostProcessDataBase(obj,mesh,fileName)
            dI.mesh    = mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            d = ps.create();
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
        
        function initPrincipalDirections(obj)
            if isempty(obj.designVariable.alpha)
                dim = obj.physicalProblem.getDimensions();
                nelem = size(obj.dvolu,1);
                ndim = dim.ndimf; %dim.ndim
                alpha0 = zeros(ndim,nelem);
                alpha0(1,:) = 1;
                obj.designVariable.alpha = alpha0;
            end
            obj.physicalProblem.variables.principalDirections = obj.designVariable.alpha;
        end
    
    end

    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        computeGradientValue(obj)
        solveState(obj)
        solveAdjoint(obj)
    end
end