classdef ShFunWithElasticPdes < Shape_Functional
    
    properties (Access = protected)
        interpolation
        matProps
        physProb
    end
    
    methods (Access = public)
     
        function computeCostAndGradient(obj,x)
            obj.updateMaterialProperties(x);
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
        
        function obj = ShFunWithElasticPdes(settings)
            obj.init(settings);
            matProp = settings.TOL;
            matType = settings.material;
            interp  = settings.method; 
            pdim    = settings.pdim;
            obj.createMaterialInterpolation(matProp,matType,interp,pdim);
        end
        
        function createEquilibriumProblem(obj,fileName)
            obj.physProb = FEM.create(fileName);
            obj.physProb.preProcess;
        end
        
        function updateMaterialProperties(obj,x)
            rho = obj.filter.getP0fromP1(x);
            obj.matProps = obj.interpolation.computeMatProp(rho);
        end
        
        function computeGradient(obj)
            g = zeros(obj.physProb.geometry.interpolation.nelem,obj.physProb.element.quadrature.ngaus);
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
        
    end
    
    methods (Access = private)
        
        function createMaterialInterpolation(obj,matProp,matType,interp,pdim)
            d.dim=pdim;
            d.interpolation=interp;
            d.typeOfMaterial=matType;
            d.constitutiveProperties=matProp;
            obj.interpolation = Material_Interpolation.create(d);
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        updateGradient(obj)
        solvePDEs(obj)
    end
end

