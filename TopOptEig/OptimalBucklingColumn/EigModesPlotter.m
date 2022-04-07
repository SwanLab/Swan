classdef EigModesPlotter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        iter
        E1
        E2
        mode1
        mode2
    end
    
    properties (Access = private)
        nElem
        length
    end
    
    methods (Access = public)
        
        function obj = EigModesPlotter(cParams)
            obj.init(cParams)            
        end
        
        function plot(obj,A,m1,m2,iter,D)
           obj.storeConvergenceVariables(iter,D)
           obj.plotColumnArea(A);
           obj.plotBucklingModes(m1,m2);
           obj.plotEigenvaluesAndIterations();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.length = cParams.length;
        end

        function storeConvergenceVariables(obj,iter,D)
            obj.iter     = iter;
            obj.E1(iter) = D(1,1);
            obj.E2(iter) = D(2,2);
        end           
        
        function plotColumnArea(obj,A)            
            L = obj.length;
            ch = 0:L:1-L;
            z = sqrt(A);
            figure(1)
            subplot(2,2,[1 3]);plot(ch,z)
            grid on
            grid minor
            title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
        end

        function plotBucklingModes(obj,m1,m2)
            L = obj.length;
            h  = 0:L:1;
            mod1 = m1;
            mod2 = m2;
            subplot(2,2,2); plot(h,-mod1(1:2:end));
            grid on
            grid minor
            title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            subplot(2,2,4); plot(h,-mod2(1:2:end));
            grid on
            grid minor
            title('Second Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
        end

        function plotEigenvaluesAndIterations(obj)
            figure(2)
            hold on
            plot(1:obj.iter,obj.E1);
            plot(1:obj.iter,obj.E2);
            hold off
            grid on
            grid minor
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Eigenvalues','Interpreter', 'latex','fontsize',18,'fontweight','b');
            axis([0 obj.iter 0 100]);
        end
        
        
    end
    
end