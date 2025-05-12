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
        velFunMat
        PFunMat
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
            plotFun = @(s, i) AirfoilOptimizer.plotAirfoilContour(s, i);
            VideoGenrator.compute('AirfoilOptimization', 1501, obj.ParamsMat, plotFun);
       end

       function generateVelVideo(obj)
            plotFun = @(s, ~) TestNaca.plotVelocity(s);
            VideoGenrator.compute('AirfoilOptimization-Velocity', 150, obj.velFunMat, plotFun);
       end

       function generatePVideo(obj)
            plotFun = @(s, ~) TestNaca.plotPressure(s);
            VideoGenrator.compute('AirfoilOptimization-Pressure', 150, obj.PFunMat, plotFun);
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
            projected(1)  = projected(1) * (projected(2) > 1e-6);
       end

       function computeOptimization(obj)
            diff     = 1;
            iter     = 1;
            maxIter  = 1e4;
        
            rho      = 0.9;            
            epsilon  = 1e-8;           
            cache    = zeros(size(obj.optimalParams));
        
            while diff > obj.tol && iter < maxIter
                gradient = obj.optimizer.computeGradient(obj.optimalParams);
        
                cache = rho * cache + (1 - rho) * (gradient.^2);
        
                adjustedGradient = obj.learningRate * gradient ./ (sqrt(cache) + epsilon);
                updatedParams    = obj.optimalParams + adjustedGradient;
        
                obj.optimalParams = obj.projectParams(updatedParams);
        
                obj.E(end + 1) = obj.optimizer.computeOutputValues(obj.optimalParams);
        
                diff = max(abs(obj.ParamsMat(end) - obj.optimalParams));
        
                obj.ParamsMat(end + 1,:) = obj.optimalParams;

                if mod(iter,10) == 0
                     [obj.velFunMat(end + 1), obj.PFunMat(end + 1)] = computeVelPFun(obj);
                end

        
                iter = iter + 1;
            end
       end

       function computeOptimization2(obj)
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

               % if mod(iter,10) == 0
               %      [obj.velFunMat(end + 1), obj.PFunMat(end + 1)] = computeVelPFun(obj);
               % end

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

       function [velFun,PFun] = computeVelPFun(obj)
            Naca.length = 8;
            Naca.height = 4;
            Naca.nx     = 420;
            Naca.M      = obj.optimalParams(1);
            Naca.p      = obj.optimalParams(2);
            Naca.t      = obj.optimalParams(3);
            Naca.chord  = 1;
            Naca.AoA    = obj.optimalParams(4);
            
            NacaClass = TestNaca(Naca);
            NacaClass.compute(); 
            velFun = NacaClass.velocityFun;
            PFun   = NacaClass.pressureFun;
            
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