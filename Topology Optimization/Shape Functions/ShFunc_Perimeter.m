classdef ShFunc_Perimeter < ShapeFunctional
    
    properties (GetAccess = public, SetAccess = private)
       filteredDensity  

    end
    
    properties (Access = protected)
        epsilon
        regularizedDensityProjection
        axes
    end
    
    properties (Access = private)
       quad
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Perimeter(cParams)
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
            obj.computePerimeterValue();
            obj.normalizeFunction();
        end
        
        function fP = addPrintableVariables(obj)
            fP{1}.value = obj.gradient;
            fP{2}.value = obj.filteredDensity.fValues;
            fP{3}.value = obj.computePerimeterIntegrandP0();
            fP{4}.value = obj.computePerimeterIntegrandP1();
        end

        function [fun, funNames] = getFunsToPlot(obj)
            mesh = obj.designVariable.mesh;

            aa.mesh    = mesh;
            aa.fValues = obj.gradient;
            gradFun = P1Function(aa);

            aa.fValues = obj.filteredDensity.fValues;
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
            f   = obj.obtainDomainFunction();
            rho = obj.filter.compute(f,'QUADRATICMASS');
            obj.filteredDensity = rho;
        end
        
        function per0 = computePerimeterIntegrandP0(obj)
            vfrac = obj.designVariable.computeVolumeFraction();
            s.mesh    = obj.designVariable.mesh;
            s.fValues = 2/(obj.epsilon)*(1 - obj.filteredDensity.fValues);
            f = P1Function(s);
            q = Quadrature.set(obj.designVariable.mesh.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            per = f.evaluate(xV);
            per = per.*vfrac;
            per0 = per;
        end
        
        function per = computePerimeterIntegrandP1(obj)
            per = 2/(obj.epsilon)*((1 - obj.filteredDensity.fValues).*(obj.Msmooth\obj.regularizedDensityProjection));
            %value2 =  sum(obj.Msmooth*per);
        end

        function computePerimeterValue(obj)
            rho        = obj.obtainDomainFunction();
            rhoe       = obj.filteredDensity;
            rhoei      = rhoe.fValues;
            s.fValues  = 1-rhoei;
            s.mesh     = obj.designVariable.mesh;
            f          = P1Function(s);
            i.type     = 'ScalarProduct';
            i.quadType = 'QUADRATICMASS';
            i.mesh     = obj.designVariable.mesh;
            int        = Integrator.create(i);
            result     = int.compute(f,rho);
            per        = 2/(obj.epsilon)*result;
            obj.value  = per;
        end
        
        function computeGradient(obj)
            obj.computeContinousGradient();
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
            pT2 = u.computePerimeter();
            pT2 - pT%
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
        
        function computeContinousGradient(obj)
            obj.gradient = 2/obj.epsilon*(1 - 2*obj.filteredDensity.fValues);
        end
        
        function computeDiscreteGradient(obj)
            %    obj.gradient = obj.Msmooth*obj.gradient;
        end
        
    end
end