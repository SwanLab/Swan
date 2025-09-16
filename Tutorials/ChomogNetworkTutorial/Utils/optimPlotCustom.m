function stop = optimPlotCustom(~,optimValues,state,varargin)
    stop = false;
    switch state
        case "init"
            grid minor
            grid on
        case "done"
            yyaxis left
            ylabel('Cost Function Value')
            yyaxis right
            ylabel('Constraint Violation')
            set(gca, 'YScale', 'log')
            ax = gca;
            ax.YAxis(1).Color = 'b';
            ax.YAxis(2).Color = 'r';
            xlabel('Iteration Count')
        otherwise
            hold on
            yyaxis right
            scatter(optimValues.funccount, optimValues.constrviolation, 20, 'red', 'filled', 'd')
            drawnow
            yyaxis left
            scatter(optimValues.funccount, optimValues.fval, 20, 'blue', 'filled', 'd')
            drawnow
            hold off
    end
end