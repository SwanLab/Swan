classdef EigModesPlotter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
      %  nElem
      %  length
       % e
       mode1
       mode2
    end
    
    methods (Access = public)
        
        function obj = EigModesPlotter(cParams)
            obj.init(cParams)            
        end
        
        function plot(obj,x,N,L,v1,v2,e,E1,E2)
           obj.computeBucklingModes(N,v1,v2);
           obj.plotEigModes(x,N,L,e,E1,E2);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
        function plotEigModes(obj,x,N,L,e,E1,E2)
            mode1 = obj.mode1;
            mode2 = obj.mode2;
        %N = obj.nElem;
            %L = obj.length;
            % axis and profile
            ch = 0:L:1-L;
            h  = 0:L:1;
            z = sqrt(x(1:N));
            % Plotting Clamped-clamped configuration
            figure(1)
            subplot(2,2,[1 3]);plot(ch,z)
            title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
            % Buckling modes
            subplot(2,2,2); plot(h,-mode1(1:2:2*N+2));
            title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            subplot(2,2,4); plot(h,-mode2(1:2:2*N+2));
            title('Second Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            figure(2)
            hold on
            plot(e,E1);
            plot(e,E2);
            hold off
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Eigenvalues','Interpreter', 'latex','fontsize',18,'fontweight','b');
            axis([0 65 0 100]);
        end
        
        function computeBucklingModes(obj,nElem,v1,v2)
            N = nElem;
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