% residual with equality constants and inequality slack variables
function [] = bp_contour(bp,x,dx)

switch bp.prob
    case(0:2)
        % don't generate contour plots
    case(2)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(0:0.1:10,0:0.1:15);
        %     min   sum (x(i)-i)^2
        %     s.t.  sum(x(i)) = 10
        
        % equation (0.1 * x1 - x2 > 1)
        eq1 = x1 + x2 - 10;
        obj = (x1-1).^2 + (x2-2).^2;
        
        bobj = obj - bp.mu * log(x1) ...
                   - bp.mu * log(x2) ...
                   - bp.mu * log(10-x1) ...
                   - bp.mu * log(100-x2);
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',250);
        hold on;
        [C,h] = contour(x1,x2,bobj);
        clabel(C,h,'Labelspacing',250);
        title('Problem 2 Contour Plot');
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq1,[0.0,0.0],'r-','LineWidth',1);
        clabel(C,h,'Labelspacing',250);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','Constraint','New Step','Start','End');       
    case(3)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(-12:0.1:5,-10:0.1:5);
        
        % equation (0.1 * x1 - x2 > 1)
        eq1 = 0.1 * x1 - x2;
        obj = x1.^2 - 2 * x1.*x2 + 4*x2.^2;
        bobj = obj - bp.mu * log(eq1-1.0) ...
                   - bp.mu * log(x1+10) ...
                   - bp.mu * log(x1+10) ...
                   - bp.mu * log(10-x1) ...
                   - bp.mu * log(10-x2);
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',250);
        hold on;
        [C,h] = contour(x1,x2,bobj);
        clabel(C,h,'Labelspacing',250);
        title('Problem 3 Contour Plot');
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq1,[1,1.1],'r-','LineWidth',1);
        clabel(C,h,'Labelspacing',250);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','Constraint','New Step','Start','End');
    case(4)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(1:0.05:5,1:0.05:5);
        
        % equation 1
        eq1 = x1 .* x2 - 5;
        % equation 2
        eq2 = -x1.^2 - x2.^2 + 20;
        % objective
        obj = x2.*(5+x1);
        for i = 1:size(eq1,1),
            for j = 1:size(eq1,2),
                if (eq1(i,j)<=0.01||eq2(i,j)<=0.01||x1(i,j)>=4.95||x1(i,j)<=1.05||x2(i,j)>=4.95||x2(i,j)<=1.05),
                    bobj(i,j) = 0;
                else
                    bobj(i,j) = obj(i,j) - bp.mu * log(eq1(i,j)) ...
                        - bp.mu * log(eq2(i,j)) ...
                        - bp.mu * log(x1(i,j)-1) - bp.mu * log(x2(i,j)-1) ...
                        - bp.mu * log(5-x1(i,j)) - bp.mu * log(5-x2(i,j));
                end
            end
        end
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',250);
        hold on;
        lines = linspace(min(min(bobj))+1,max(max(bobj))-1,20);
        [C,h] = contour(x1,x2,bobj,lines);
        clabel(C,h,'Labelspacing',250);
        title('Problem 4 Contour Plot');
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq1,[0,0.1],'r-','LineWidth',1);
        clabel(C,h,'Labelspacing',250);
        [C,h] = contour(x1,x2,eq2,[0,0.1],'b-','LineWidth',1);
        clabel(C,h,'Labelspacing',250);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','Constraint','New Step','Start','End');
    case(5)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(-12:0.1:5,-10:0.1:5);
        
        % equation (0.1 * x1 - x2 > 1)
        eq1 = 0.1 * x1 - x2;
        % equation (-10*x1 + x2 > 1)
        eq2 = -10 * x1 + x2;
        obj = x1.^2 - 2 * x1.*x2 + 4*x2.^2;
        for i = 1:size(eq1,1),
            for j = 1:size(eq1,2),
                if (eq1(i,j)<=1.01||eq2(i,j)<=1.01),
                    bobj(i,j) = 0;
                else
                    bobj(i,j) = obj(i,j) - bp.mu * log(eq1(i,j)-1.0) ...
                        - bp.mu * log(eq2(i,j)-1.0);
                end
            end
        end
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',250);
        hold on;
        lines = linspace(min(min(bobj))+1,max(max(bobj))-1,50);
        [C,h] = contour(x1,x2,bobj,lines);
        clabel(C,h,'Labelspacing',250);
        title('Problem 5 Contour Plot');
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq1,[1,1.1],'r-','LineWidth',1);
        clabel(C,h,'Labelspacing',250);
        [C,h] = contour(x1,x2,eq2,[1,1.1],'b-','LineWidth',1);
        clabel(C,h,'Labelspacing',250);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','Constraint','New Step','Start','End');
    case(6)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(0:0.05:5,0:0.05:10);
        
        % equation (2 * x1 + x2 <= 9)
        eq1 = 2 * x1 + x2;
        % equation (x1 + 2 * x2  = 10)
        eq2 = x1 + 2 * x2;
        obj = x1.^2 + 2 * x2.^2;
        for i = 1:size(eq1,1),
            for j = 1:size(eq1,2),
                if ((9-eq1(i,j))<=0.01 || x1(i,j)<0.01 || x2(i,j)<0.01),
                    bobj(i,j) = 0;
                else
                    bobj(i,j) = obj(i,j) - bp.mu * log(9 - eq1(i,j)) ...
                        - bp.mu * log(x1(i,j)) - bp.mu * log(x2(i,j));
                end
            end
        end
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',500);
        hold on;
        lines = linspace(min(min(bobj))+1,max(max(bobj))-1,30);
        [C,h] = contour(x1,x2,bobj,lines);
        clabel(C,h,'Labelspacing',500);
        title('Problem 6 Contour Plot');
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq1,[9,9],'r--','LineWidth',2);
        clabel(C,h,'Labelspacing',500);
        [C,h] = contour(x1,x2,eq2,[10,10],'g-.','LineWidth',2);
        clabel(C,h,'Labelspacing',500);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','2*x_1+x_2<=9','x1+2*x2=10','New Step','Start','End');
    case(7)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(0:0.1:10,0:0.1:10);
        
        % equation (x1 + x2  = 10)
        eq = x1 + x2;
        obj = (x1-5).^2 + (x2-5).^2;
        for i = 1:size(eq,1),
            for j = 1:size(eq,2),
                if (x1(i,j)<0.01 || x2(i,j)<0.01 || x1(i,j)>9.99 || x2(i,j)>8.99),
                    bobj(i,j) = 0;
                else
                    bobj(i,j) = obj(i,j) ...
                        - bp.mu * log(x1(i,j)) - bp.mu * log(x2(i,j)) ...
                        - bp.mu * log(10 - x1(i,j)) - bp.mu * log(9 - x2(i,j));
                end
            end
        end
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',500);
        hold on;
        lines = linspace(min(min(bobj))+1,max(max(bobj))-1,30);
        [C,h] = contour(x1,x2,bobj,lines);
        clabel(C,h,'Labelspacing',500);
        title(['Problem 7 Contour Plot with mu = ' num2str(bp.mu)]);
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq,[10,10],'r-.','LineWidth',2);
        clabel(C,h,'Labelspacing',500);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','x1+x2=10','New Step','Start','End');
    case(8)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(0:0.1:10,0:0.1:10);
        
        % equation x1^2 <= 9 converted to 9 - x1^2 - x2
        eq = 9 - x1.^2 - x2;
        obj = (x1-5).^2;
        for i = 1:size(eq,1),
            for j = 1:size(eq,2),
                if (x1(i,j)<0.01 || x2(i,j)<0.01 || x1(i,j)>9.99),
                    bobj(i,j) = 0;
                else
                    bobj(i,j) = obj(i,j) ...
                        - bp.mu * log(x1(i,j)) - bp.mu * log(x2(i,j)) ...
                        - bp.mu * log(10 - x1(i,j));
                end
            end
        end
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',500);
        hold on;
        lines = linspace(min(min(bobj))+1,max(max(bobj))-1,30);
        [C,h] = contour(x1,x2,bobj,lines);
        clabel(C,h,'Labelspacing',500);
        title(['Problem 8 Contour Plot with mu = ' num2str(bp.mu)]);
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq,[0,0],'r-.','LineWidth',2);
        clabel(C,h,'Labelspacing',500);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','9 - x1^2 - x2=10','New Step','Start','End');
    case(10)
        %% generate contour plot
        % design variables at mesh points
        [x1,x2] = meshgrid(-8:0.1:5,-8:0.1:5);
        
        % equation (0.1 * x1 - x2 > 1)
        eq1 = 0.1 * x1 - x2;
        obj = x1.^2 - 2 * x1.*x2 + 4*x2.^2;
        for i = 1:size(eq1,1),
            for j = 1:size(eq1,2),
                if (eq1(i,j)<=1.01),
                    bobj(i,j) = 0;
                else
                    bobj(i,j) = obj(i,j) - bp.mu * log(eq1(i,j)-1.0);
                end
            end
        end
        
        figure(100)
        hold off;
        [C,h] = contour(x1,x2,obj);
        clabel(C,h,'Labelspacing',250);
        hold on;
        lines = linspace(min(min(bobj))+1,max(max(bobj))-1,50);
        [C,h] = contour(x1,x2,bobj,lines);
        clabel(C,h,'Labelspacing',250);
        title('Problem 10 Contour Plot');
        xlabel('x_1');
        ylabel('x_2');
        hold on;
        % solid lines to show constraint boundaries
        [C,h] = contour(x1,x2,eq1,[1,1.1],'r-','LineWidth',1);
        clabel(C,h,'Labelspacing',250);
        % plot step
        plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
        plot(x(1),x(2),'ro');
        plot(x(1)+dx(1),x(2)+dx(2),'bx');
        % show a legend
        legend('Objective','Barrier Problem','Constraint','New Step','Start','End');
end

