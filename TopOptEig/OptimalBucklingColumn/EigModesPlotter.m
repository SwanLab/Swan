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
        
        function plot(obj,x,v1,v2,iter,D)
           obj.computeConvergence2Eigenvalues(iter,D)
           obj.computeBucklingModes(v1,v2);
           obj.plotColumnArea(x);
           obj.plotBucklingModes();
           obj.plotEigenvaluesAndIterations();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem  = cParams.nElem;
            obj.length = cParams.length;
        end

        function computeConvergence2Eigenvalues(obj,iter,D)
            obj.iter     = iter;
            obj.E1(iter) = D(1,1);
            obj.E2(iter) = D(2,2);
        end           
        
        function plotColumnArea(obj,x)
            N = obj.nElem;
            L = obj.length;
            ch = 0:L:1-L;
            z = sqrt(x(1:N));
            figure(1)
            subplot(2,2,[1 3]);plot(ch,z)
            grid on
            grid minor
            title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
        end

        function plotBucklingModes(obj)
            N = obj.nElem;
            L = obj.length;
            h  = 0:L:1;
            mod1 = obj.mode1;
            mod2 = obj.mode2;
            subplot(2,2,2); plot(h,-mod1(1:2:2*N+2));
            grid on
            grid minor
            title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            subplot(2,2,4); plot(h,-mod2(1:2:2*N+2));
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
            axis([0 65 0 100]);
        end
        
        function computeBucklingModes(obj,v1,v2)
            N = obj.nElem;
            Mode1=zeros(2*N+2);
            Mode2=zeros(2*N+2);
            for i=3:2*N
                Mode1(i)=v1(i-2);
                Mode2(i)=v2(i-2);
            end
            obj.mode1 = Mode1;
            obj.mode2 = Mode2;
        end         
        
        
    end
    
end