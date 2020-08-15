classdef ShFunc_Perimeter < ShapeFunctional
    
    properties (GetAccess = public, SetAccess = private)
       regularizedDensity        
    end
    
    properties (Access = protected)
        epsilon
        regularizedDensityProjection
        axes
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Perimeter(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            cParams.filterParams.filterType = 'PDE';
            cParams.filterParams.femSettings.bcApplierType = 'Neumann';         
            obj.init(cParams);
          %  obj.initFrame();
        end
        
        function computeCostAndGradient(obj)
            obj.updateProtectedVariables()
            obj.computeRegularizedDensity()
            obj.computeRegularizedDensityProjection()
            obj.computePerimeterValue()
            obj.computePerimeterGradient()
        end
        
        function d = addPrintableVariables(obj,d)
            d.gradient   = obj.gradient;
            d.regDensity = obj.regularizedDensity;
        end        
        
    end
    
    methods (Access = protected)
        
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
             obj.regularizedDensity = obj.filter.getP1fromP1(obj.designVariable);
%              rho = obj.regularizedDensity;
%              m = obj.filter.diffReacProb.mesh;
%              
%              
%              x = m.coord(:,1);
%              y = m.coord(:,2);
%              z = rho;
%              tri = delaunay(x,y);
%              f = figure();
%              trisurf(tri,x,y,z);
%              shading interp
%              
%              node1 = rho(m.connec(:,1));
%              node2 = rho(m.connec(:,2));
%              node3 = rho(m.connec(:,3));
%              rhoElem = mean([node1,node2,node3],2);
%            cla(obj.axes)
%             patchHandle = patch('Faces',m.connec,'Vertices',m.coord,...
%                 'FaceAlpha','flat','EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
%            set(patchHandle,'FaceVertexAlphaData',rhoElem,'FaceAlpha','flat');
          end
        
        function initFrame(obj)
            figHandle = figure();
            
            set(figHandle,'Pointer','arrow','NumberTitle','off');
            
            hold on
            axis off
            axis equal
            
            obj.axes = figHandle.Children;
        end
        
        function computeRegularizedDensityProjection(obj)
            obj.regularizedDensityProjection = obj.filter.integrate_L2_function_with_shape_function(obj.designVariable);
        end
        
        function computePerimeterValue(obj)
            obj.value = 2/(obj.epsilon)*((1 - obj.regularizedDensity)'*obj.regularizedDensityProjection);
        end
        
        function computePerimeterGradient(obj)
            obj.computeContinousGradient();
            obj.computeDiscreteGradient();
        end
        
        function computeContinousGradient(obj)
            obj.gradient = 2/obj.epsilon*(1 - 2*obj.regularizedDensity);
        end
        
        function computeDiscreteGradient(obj)
          %  obj.gradient = obj.Msmooth*obj.gradient;
        end
        
    end
end