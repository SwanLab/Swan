classdef ShFunWithElasticPdes < ShapeFunctional
    
    properties (Access = protected)
        interpolation
        physicalProblem
    end
    
    methods (Access = public)
        
        function computeCostAndGradient(obj)
            obj.nVariables = obj.designVariable.nVariables;
            obj.updateHomogenizedMaterialProperties();
            obj.solvePDEs();
            obj.computeFunctionValue();
            obj.computeGradient();
            obj.normalizeFunctionAndGradient()
        end
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
        end
        
    end
    
    methods (Access = protected)
        
        function createEquilibriumProblem(obj,fileName)
            obj.physicalProblem = FEM.create(fileName);
            obj.physicalProblem.preProcess;
            obj.initPrincipalDirections();
        end
        
        function updateHomogenizedMaterialProperties(obj)
            nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
            x  = obj.designVariable.value;
            for ivar = 1:obj.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xf(:,ivar) = obj.filter.getP0fromP1(xs);
            end                        
            alpha = obj.physicalProblem.variables.principalDirections;
            obj.homogenizedVariablesComputer.computeCtensor(xf,alpha);
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
            ndim = 2;
            nelem = obj.physicalProblem.element.nelem;
            obj.physicalProblem.variables.principalDirections = zeros(ndim,ndim,nelem);
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        updateGradient(obj)
        solvePDEs(obj)
    end
end