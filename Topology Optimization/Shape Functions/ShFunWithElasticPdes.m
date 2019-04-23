classdef ShFunWithElasticPdes < Shape_Functional
    
    properties (Access = protected)
        interpolation
        homogenizedMatProps
        physProb
        homogenizedVariablesComputer
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
            f{1} = obj.physProb;
        end
        
    end
    
    methods (Access = protected)
        
        function createEquilibriumProblem(obj,fileName)
            obj.physProb = FEM.create(fileName);
            obj.physProb.preProcess;
        end
        
        function updateHomogenizedMaterialProperties(obj,x)
            rho = obj.filter.getP0fromP1(x);
            obj.homogenizedVariablesComputer.computeMatProp(rho);
        end
        
        function computeGradient(obj)
            nelem = obj.elemGradientSize.nelem;
            ngaus = obj.elemGradientSize.ngaus;
            g = zeros(nelem,ngaus);
            for igaus = 1:obj.physProb.element.quadrature.ngaus
                for istre = 1:obj.physProb.element.getNstre()
                    for jstre = 1:obj.physProb.element.getNstre()
                        g(:,igaus) = g(:,igaus) + obj.updateGradient(igaus,istre,jstre);
                    end
                end
            end
            g = obj.filter.getP1fromP0(g);
            g = obj.Msmooth*g;
            obj.gradient = g;
        end
        
        function createHomogenizedVariablesComputer(obj,cParams)
            cP = cParams.materialInterpolationParams;
            cP.nelem = obj.elemGradientSize.nelem;
            cP.ngaus = obj.elemGradientSize.ngaus;
            h = HomogenizedVarComputer.create(cP);
            obj.homogenizedVariablesComputer = h;
        end
        
    end

    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        updateGradient(obj)
        solvePDEs(obj)
    end
end

