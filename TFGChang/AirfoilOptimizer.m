classdef AirfoilOptimizer < handle

    properties (Access = public)
        optimalParams
    end
    
    properties (Access = private)
        features
        learningRate
        optimizer
        tol
        upperBC
        lowerBC
        ParamsMat
        E
    end
  
   methods (Access = public)

       function obj = AirfoilOptimizer(cParams)
           obj.init(cParams);
       end

       function computeOptAirfoilParams(obj)
           obj.computeOptimization();
       end

       function plotEEvolution(obj)
           obj.plotE();
       end

       function generateAFSOPVideo(obj)
           obj.generateVideo();
       end

       function saveData(obj)
            save("OptData.mat","obj");
       end


   end

   methods (Access = private)

       function init(obj,cParams)
           obj.features         = cParams.features;
           obj.optimizer        = cParams.optimizer;
           obj.tol              = cParams.tol;
           obj.learningRate     = cParams.learningRate;
           obj.upperBC          = cParams.upperBC;
           obj.lowerBC          = cParams.lowerBC;
           obj.optimalParams    = obj.features; 
           obj.E(1)             = obj.optimizer.computeOutputValues(obj.optimalParams);
           obj.ParamsMat(1,:)   = obj.features;
       end  

       function projected = projectParams(obj,params)
            projected     = max(min(params, obj.upperBC), obj.lowerBC);
            projected(2)  = projected(2) * (projected(1) > 1e-6);
       end

       function computeOptimization(obj)
           diff     = 1;
           iter     = 1;
           maxIter  = 1e4;

           while diff > obj.tol && iter < maxIter

               % if (iter > 130)
               %     obj.learningRate = 1;
               % end

               gradient          = obj.optimizer.computeGradient(obj.optimalParams);
               updatedParams     = obj.optimalParams + obj.learningRate * gradient;
               obj.optimalParams = obj.projectParams(updatedParams);

               obj.E(end + 1)      = obj.optimizer.computeOutputValues(obj.optimalParams);
                             
               diff                = max(abs(obj.ParamsMat(end) - obj.optimalParams));
               %max(abs(obj.ParamsMat(end) - obj.optimalParams));
               %max(abs(gradient));
               %max(abs(obj.E(end - 1) - obj.E(end)));
               obj.ParamsMat(end + 1,:) = obj.optimalParams;

               iter = iter + 1;
           end

       end

       function plotE(obj)
            figure;
            plot(1:length(obj.E), obj.E);
            xlabel('Iteration');
            ylabel('Aerodynamic Efficiency');
            title('Evolution of Airfoil Aerodynamic Efficiency vs Optimization Iteration');
            grid on;
       end

       function generateVideo(obj)

            v = VideoWriter('AirfoilOptimization');
            v.FrameRate = 30;
            open(v);
     
            figure;
            for i = 1:1501
                s = obj.ParamsMat(i,:);
                AirfoilOptimizer.plotAirfoilContour(s,i)
                frame = getframe(gcf);
    
                % Controlar la velocidad según el frame actual
                if i <= 130
                    % Desacelerar: mostrar cada frame 3 veces
                    repeat = 2;
                else
                    % Acelerar: saltar algunos frames (solo guardar uno de cada 3)
                    if mod(i,6) ~= 0
                        continue
                    end
                    repeat = 1;
                end
                
                for r = 1:repeat
                    writeVideo(v, frame);
                end
            end
            
            close(v);
        end
         
   end

   methods (Static)

        function [xContour,yContour] = computeAirfoilContour(s)
            x   = 0:0.001:1;
            m   = s(1);
            p   = s(2);
            t   = s(3);

            if (m > 1e-6)
                yc   = (x<=p).*(m/p^2.*(2*p*x - x.^2)) + ...
                       (x>p).*(m/(1-p)^2.*((1 - 2*p) + 2*p*x - x.^2));
                dydx = (x<=p).*(2*m/p^2.*(p - x)) + ...
                       (x>p).*(2*m/(1-p)^2.*(p - x));
            else
                yc   = 0;
                dydx = 0;
            end

            yt   = 5*t.*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1036*x.^4);
            theta = atan(dydx);
            
            xu = x - yt.*sin(theta);
            yu = yc + yt.*cos(theta);
            xl = x + yt.*sin(theta);
            yl = yc - yt.*cos(theta);
            
            xContour = [xu, fliplr(xl)];
            yContour = [yu, fliplr(yl)];
       end

       function [xRot,yRot] = rotateAirfoil(xContour,yContour,AoA)
            xRot = xContour*cos(AoA) - yContour*sin(AoA);
            yRot = xContour*sin(AoA) + yContour*cos(AoA);
       end

       function plotAirfoilContour(s,i)
            
            AoA = -deg2rad(s(4));

            [xContour,yContour] = AirfoilOptimizer.computeAirfoilContour(s);
            
            [xRot,yRot] = AirfoilOptimizer.rotateAirfoil(xContour,yContour,AoA);
            
            plot(xRot, yRot, 'k-', 'LineWidth', 1);
            xlim([-0.2 1.2]);  
            ylim([-0.6 0.6]); 
            axis equal; 
            xlabel('x'); ylabel('y');
            title(sprintf('Airfoil Shape Optimization - Iteration %d', i - 1));
            grid on;

       end


   end

end