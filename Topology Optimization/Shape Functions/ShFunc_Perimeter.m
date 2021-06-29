classdef ShFunc_Perimeter < ShapeFunctional
    
    properties (GetAccess = public, SetAccess = private)
       regularizedDensity  

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
            cParams.filterParams.quadratureOrder = 'LINEAR';
            cParams.filterParams.filterType = 'PDE';
            cParams.filterParams.femSettings.bcApplierType = 'Neumann';         
            obj.init(cParams);
          %  obj.initFrame();
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
            s.connec = obj.designVariable.mesh.connec;
            s.type   = obj.designVariable.mesh.type;
            s.fNodes = 2/(obj.epsilon)*(1 - obj.regularizedDensity);
            f = FeFunction(s);
            per = f.computeValueInCenterElement();
            per = per.*vfrac;
            per0 = per;
        end        
        
        function per = computePerimeterIntegrandP1(obj)
            per = 2/(obj.epsilon)*((1 - obj.regularizedDensity).*(obj.Msmooth\obj.regularizedDensityProjection));
            %value2 =  sum(obj.Msmooth*per);
        end
        
        function computePerimeterValue(obj)
            int = 2/(obj.epsilon)*((1 - obj.regularizedDensity).*obj.regularizedDensityProjection);
            obj.value =  sum(int);    
        end
        
        function computeGradient(obj)
            obj.computeContinousGradient();
            obj.computeDiscreteGradient();
            obj.normalizeGradient();
        end        
        
        function pT = computeGeometricTotalPerimeter(obj)
            if isequal(class(obj.designVariable),'LevelSet')            
            u = obj.designVariable.getUnfittedMesh();
            pR = obj.computeGeometricRelativePerimeter();
            pB = u.unfittedBoundaryMesh.computeVolume();            
            pT = pR + pB;
            pT2 = u.computePerimeter();
            pT2 - pT
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
                p = u.boundaryCutMesh.mesh.computeVolume;
            else
                p = 0;
            end
        end
        
        function computeContinousGradient(obj)
            obj.gradient = 2/obj.epsilon*(1 - 2*obj.regularizedDensity);
        end
        
        function computeDiscreteGradient(obj)
            %    obj.gradient = obj.Msmooth*obj.gradient;
        end
        
        
        
    end
end