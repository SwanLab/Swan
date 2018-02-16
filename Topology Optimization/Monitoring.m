classdef Monitoring < handle
    %Monitoring Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        figures
        monitor
        nfigs
        ncost
        nconstraint
        nstop
        stop_names
    end
    
    methods
        function obj = Monitoring(settings)
            obj.getStopVarsNames(settings.optimizer);
            
            obj.figures = {};
            obj.createFigure('Cost');
            
            obj.ncost = length(settings.cost);
            for i = 1:obj.ncost
                if isempty(settings.multipliers)
                    obj.createFigure([obj.setCase(settings.cost{i}) ' (wt. 1.0)']);
                else
                    obj.createFigure([obj.setCase(settings.cost{i}) sprintf(' (wt. %.2f)',settings.multipliers(i))]);
                end
            end
            
            obj.nconstraint = length(settings.constraint);
            for i = 1:obj.nconstraint
                obj.createFigure(['Cstr ' num2str(i) ': ' obj.setCase(settings.constraint{i})]);
                obj.createFigure(['\lambda' obj.subindex(obj.setCase(settings.constraint{i}))]);
            end
            
            for i = 1:obj.nstop
                obj.createFigure(['Conv Criteria ' num2str(i) ': ' obj.stop_names{i}]);
            end
            
            obj.nfigs = length(obj.figures);
            if obj.nfigs <= 4
                nrows = 1;
            else
                nrows = 2;
            end
            ncols = round(obj.nfigs/nrows);
            
            % Create obj.figures
            obj.monitor = figure; hold on
            set(obj.monitor,'units','normalized','outerpos',[0 0 1 1]);
            for i = 1:obj.nfigs
                obj.figures{i}.style = subplot(nrows,ncols,i);
                switch obj.figures{i}.chart_type
                    case 'plot'
                        obj.figures{i}.handle = plot(0,0);
                        title(obj.figures{i}.title);
                        grid on
                    case 'bar'
                        obj.figures{i}.handle = bar(0,0);
                        title(obj.figures{i}.title);
                        grid on
                end
            end
        end
        
        function display(obj,iteration,cost,constraint,lambda,stop_vars)
            set(obj.monitor,'NumberTitle','off','Name',sprintf('Monitoring - Iteration: %.0f',iteration))
            obj.figures{1} = obj.updateFigure(obj.figures{1},iteration,cost.value);
            for i = 1:obj.ncost
                k = i+1;
                obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,cost.ShapeFuncs{i}.value);
            end
            for i = 1:obj.nconstraint
                k = i*2+obj.ncost;
                obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,constraint.ShapeFuncs{i}.value);
                k = k+1;
                obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,lambda(i));
            end
            for i = 1:obj.nstop
                k = i+1+obj.ncost+2*obj.nconstraint;
                obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,stop_vars(i,1));
            end
        end
    end
    methods (Access = private)
        function obj = createFigure(obj,TITLE)
            obj.figures{end+1}.title = TITLE;
            obj.figures{end}.iteration = [];
            obj.figures{end}.variable = [];
            if contains(TITLE,'kappa')
                obj.figures{end}.chart_type = 'bar';
            else
                obj.figures{end}.chart_type = 'plot';
            end
        end
        
        function obj = getStopVarsNames(obj,optimizer)
            switch optimizer
                case 'SLERP'
                    obj.stop_names = {'\Deltacost';'\Deltavolume';'\kappa'};
                case 'PROJECTED GRADIENT'
                    obj.stop_names = {'\Deltacost';'\Deltavolume';'\kappa'};
                case 'MMA'
                    obj.stop_names = {'kktnorm';'outit'};
                case 'IPOPT'
                    obj.stop_names = {};
            end
            obj.nstop = length(obj.stop_names);
        end
    end
    methods (Access = private, Static)
        function figures = updateFigure(figures,iteration,variable)
            figures.iteration(end+1) = iteration;
            figures.variable(end+1) = variable;
            set(figures.handle,'XData',figures.iteration,'YData',figures.variable);
            set(figures.style,'XLim',[0 figures.iteration(end)])
            drawnow
        end
        
        function str = setCase(str)
            str = [upper(str(1)) lower(str(2:end))];
        end
        
        function str = subindex(str_)
            str = [];
            for i = 1:length(str_)
                str(end+1:end+2) = ['_' str_(i)];
            end
        end
    end
    
end

