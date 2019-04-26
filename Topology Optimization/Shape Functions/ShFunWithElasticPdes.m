classdef ShFunWithElasticPdes < ShapeFunctional
    
    properties (Access = protected)
        interpolation
        physicalProblem
    end
    
    methods (Access = public)
        
        function computeCostAndGradient(obj,x)
            obj.updateHomogenizedMaterialProperties(x);
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
        
        function updateHomogenizedMaterialProperties(obj,x)
            rho = obj.filter.getP0fromP1(x);
            alpha = obj.physicalProblem.variables.principalDirections;
            obj.homogenizedVariablesComputer.computeCtensor(rho,alpha);
        end
        
        function computeGradient(obj)
            nelem = obj.physicalProblem.geometry.interpolation.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            nstre = obj.physicalProblem.element.getNstre();
            g = zeros(nelem,ngaus);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        g(:,igaus) = g(:,igaus) + obj.updateGradient(igaus,istre,jstre);
                    end
                end
            end
            g = obj.filter.getP1fromP0(g);
            g = obj.Msmooth*g;
            obj.gradient = g;
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