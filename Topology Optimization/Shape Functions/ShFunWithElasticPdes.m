classdef ShFunWithElasticPdes < ShapeFunctional
    
    properties (Access = protected)
        interpolation
        physicalProblem
    end
    
    properties (Access = private)
        orientationUpdater        
    end
    
    properties (Access = protected)
        regDesignVariable
    end
    
    methods (Access = public)
        
        function computeCostAndGradient(obj)
            obj.nVariables = obj.designVariable.nVariables;
            for i = 1:1
                obj.updateHomogenizedMaterialProperties();
                obj.solvePDEs();
                obj.updateAlpha();
            end
            obj.computeFunctionValue();
            obj.computeGradient();
            obj.normalizeFunctionAndGradient()
        end
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
        end
        
        function f = getRegularizedDensity(obj)
            f = obj.regDesignVariable;
        end
        
        function d = addPrintableVariables(obj,d)
            d.phyProblems = obj.getPhysicalProblems();
            d.regDensity  = obj.getRegularizedDensity();
        end        
        
    end
    
    methods (Access = protected)
        
        function createOrientationUpdater(obj)
            cParams.type = 'MinimumEigenValue';
            obj.orientationUpdater = OrientationUpdater.create(cParams);            
        end
        
        function createEquilibriumProblem(obj,fileName)
            obj.physicalProblem = FEM.create(fileName);
            obj.initPrincipalDirections();
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
        end
        
        function filterDesignVariable(obj)
            nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
            x  = obj.designVariable.value;
            for ivar = 1:obj.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xf(:,ivar) = obj.filter.getP0fromP1(xs);
            end
            obj.regDesignVariable = xf;
        end
        
        function computeGradient(obj)
            nelem = obj.physicalProblem.geometry.interpolation.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            nstre = obj.physicalProblem.element.getNstre();
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        g(:,igaus,:) = g(:,igaus,:) + obj.updateGradient(igaus,istre,jstre);
                    end
                end
            end
            
            gf = zeros(size(obj.Msmooth,1),obj.nVariables);
            for ivar = 1:obj.nVariables
                gs = squeeze(g(:,:,ivar));
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            g = obj.Msmooth*gf;
            obj.gradient = g(:);
        end
        
    end
    
    methods (Access = private)
        
        function initPrincipalDirections(obj)
            ndim = obj.physicalProblem.mesh.ndim;
            nelem = obj.physicalProblem.element.nelem;
            alpha0 = zeros(ndim,nelem);
            alpha0(1,:) = 1;
            obj.physicalProblem.variables.principalDirections = alpha0;
            obj.designVariable.alpha = alpha0;
        end
        
        function updateAlpha(obj)
            cParams.pD = obj.physicalProblem.variables.principalDirections;
            cParams.pS = obj.physicalProblem.variables.principalStress;
            obj.orientationUpdater.compute(cParams);
            alpha = obj.orientationUpdater.alpha;            
            obj.designVariable.alpha = alpha;
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        updateGradient(obj)
        solvePDEs(obj)
    end
end