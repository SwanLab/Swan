classdef ShFunc_PerimeterWithPNorm < ShapeFunctional
    
    properties (GetAccess = public, SetAccess = private)
       regularizedDensity  
       normExponent = 8
    end
    
    properties (Access = protected)
        epsilon
        regularizedDensityProjection
        axes
    end
    
    properties (Access = private)
       quad
       valueNonDimensioned
    end
    
    methods (Access = public)
        
        function obj = ShFunc_PerimeterWithPNorm(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            cParams.filterParams.filterType = 'PDE';
            obj.init(cParams);
          obj.computeFunctionAndGradient();
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end
        
        function computeFunction(obj)
            obj.updateProtectedVariables();
            obj.computeRegularizedDensity();
            obj.computeRegularizedDensityProjection();
            obj.computePerimeterValue();
            obj.normalizeFunction();
        end
        
        function fP = addPrintableVariables(obj)
            fP{1}.value = obj.gradient;
            fP{2}.value = obj.regularizedDensity;
            fP{3}.value = obj.computePerimeterIntegrandP0();
            fP{4}.value = obj.computePerimeterIntegrandP1();
        end

        function [fun, funNames] = getFunsToPlot(obj)
            mesh = obj.designVariable.mesh.meshes{1};

            aa.mesh    = mesh;
            aa.fValues = obj.gradient;
            gradFun = P1Function(aa);

            aa.fValues = obj.regularizedDensity;
            regDensFun = P1Function(aa);

            aa.fValues = squeeze(obj.computePerimeterIntegrandP0());
            perIntegP0 = P0Function(aa);

            aa.fValues = obj.computePerimeterIntegrandP1();
            perIntegP1 = P1Function(aa);

            fun = {gradFun, regDensFun, perIntegP0, perIntegP1};
            funNames = {'PerimeterGradient', 'RegularizedDensity', ...
                        'PerimeterP0', 'PerimeterP1'};
        end
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
            v{2} = obj.computeGeometricRelativePerimeter();
            v{3} = obj.computeGeometricTotalPerimeter();
            v{4} = obj.computeEpsilonOverhValue();
        end
        
        function t = getTitlesToPlot(obj)
            t{1} = 'Perimeter non scaled';
            t{2} = 'Geometric Relative Perimeter';
            t{3} = 'Geometric Total Perimeter';
            t{4} = 'epsilon over h';
        end
        
        function fP = createPrintVariables(obj)
            types = {'ScalarNodal','ScalarNodal','ScalarGauss','ScalarNodal'};
            names = {'PerimeterGradient','RegularizedDensity',...
                     'PerimeterGauss','PerimeterNodal'};
            fP = obj.obtainPrintVariables(types,names);
        end
        
    end
    
    methods (Access = private)
        
        function updateProtectedVariables(obj)
            obj.updateEpsilonValue()
            obj.updateEpsilonInFilter()
        end
        
        function updateEpsilonValue(obj)
            obj.epsilon = obj.target_parameters.epsilon_perimeter;
        end
        
        function updateEpsilonInFilter(obj)
            obj.filter.updateEpsilon(obj.epsilon);
        end
        
        function computeRegularizedDensity(obj)
             obj.regularizedDensity = obj.filter.getP1fromP1(obj.designVariable.value);
         end
        
        function computeRegularizedDensityProjection(obj)
            obj.regularizedDensityProjection = obj.filter.integrate_L2_function_with_shape_function(obj.designVariable.value);
        end
        
        function per0 = computePerimeterIntegrandP0(obj)
            vfrac = obj.designVariable.computeVolumeFraction();
            s.mesh    = obj.designVariable.mesh;
            s.fValues = 2/(obj.epsilon)*(1 - obj.regularizedDensity);
            f = P1Function(s);
            q = Quadrature.set(obj.designVariable.mesh.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            per = f.evaluate(xV);
            per = per.*vfrac;
            per0 = per;
        end
        
        function per = computePerimeterIntegrandP1(obj)
            per = 2/(obj.epsilon)*((1 - obj.regularizedDensity).*(obj.Msmooth\obj.regularizedDensityProjection));
        end
        
        function computePerimeterValue(obj)
            p1.fValues = obj.designVariable.value;
            p1.mesh    = obj.designVariable.mesh;
            s.fun      = P1Function(p1);
            s.mesh     = obj.designVariable.mesh;
            s.exponent = obj.normExponent;
            s.type     = 'ShapeFunctionWithSomeFunctionPowered';
            s.designVarType = obj.designVariable.type;
            rhs        = RHSintegrator.create(s);

            designVarProj = rhs.compute();
            rhoe          = obj.regularizedDensity';
            p             = obj.normExponent;
            obj.value     = 2/(obj.epsilon)*((1-rhoe).^p*designVarProj)^(1/p);

            obj.valueNonDimensioned = obj.value;

            %
            ss.mesh = obj.designVariable.mesh;
            ss.fValues = ((obj.designVariable.value).*(1-obj.regularizedDensity)).^p;
            ss.filename = 'integrand';
            integrand = P1Function(ss);
            integrand.print(ss);
            %
        end

        function computeGradient(obj)
            obj.computeDiscreteGradient();
            obj.normalizeGradient();
        end        
        
        function pT = computeGeometricTotalPerimeter(obj)
            if isequal(class(obj.designVariable),'LevelSet')
            u = obj.designVariable.getUnfittedMesh();
            if ~isempty(u.boundaryCutMesh)
            pR = obj.computeGeometricRelativePerimeter();
            pB = u.unfittedBoundaryMesh.computeVolume();
            pT = pR + pB;
            else
               pT = 0;
            end
            else 
                pT = 0;
            end
        end
        
        function eh = computeEpsilonOverhValue(obj)
            e = obj.epsilon;
            m = obj.designVariable.mesh;
            h = m.computeMeanCellSize();
            eh = e/h;
        end
        
        function p = computeGeometricRelativePerimeter(obj)
            if isequal(class(obj.designVariable),'LevelSet')
                u = obj.designVariable.getUnfittedMesh();
                if ~isempty(u.boundaryCutMesh)
                p = u.boundaryCutMesh.mesh.computeVolume;
                else 
                p = 0;
                end
            else
                p = 0;
            end
        end

        function computeDiscreteGradient(obj)
            p1.fValues = obj.designVariable.value;
            p1.mesh    = obj.designVariable.mesh;
            s.fun      = P1Function(p1);
            s.mesh     = obj.designVariable.mesh;
            s.exponent = obj.normExponent;
            s.type     = 'ShapeFunctionsWithSomeFunctionPowered';
            s.designVarType = obj.designVariable.type;
            lhs        = LHSintegrator.create(s);

            designVarProjMass = lhs.compute();
            rhoe          = obj.regularizedDensity';
            p             = obj.normExponent;

            Per       = obj.valueNonDimensioned;
            factorFun = (1/Per*(1-rhoe)*2/obj.epsilon).^(p-1);
            factorDer = 2/(obj.epsilon)*(1-2*rhoe);

            DxPer        = (factorDer.*factorFun)*designVarProjMass;
            obj.gradient = DxPer';
        end

    end
end