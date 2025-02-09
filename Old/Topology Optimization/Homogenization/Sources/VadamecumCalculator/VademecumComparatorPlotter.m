classdef VademecumComparatorPlotter < handle
    
    properties (Access = private)
        comparator
        path
        indeces
        fig
        istre
        jstre
        cleg
        h
    end
    
    
    methods (Access = public)
        
        function obj = VademecumComparatorPlotter(cParams)
            obj.init(cParams)
        end
        
        function plot(obj)
           obj.plotConstitutiveTensors();
        end
        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.comparator = cParams.vadSAcomparator;
            obj.path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';            
        end
        
        function plotConstitutiveTensors(obj)
            obj.indeces = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
            for index = 1:length(obj.indeces)
                obj.fig = figure(index);
                obj.istre = obj.indeces(index,1);
                obj.jstre = obj.indeces(index,2);
                obj.plotVademecumTensor();
                obj.plotSimpAllTensor();
                obj.addLegend();
                obj.printPlot();
            end
        end
        
        function plotSimpAllTensor(obj)
            C = obj.comparator.CtensorSIMPALL;
            rho = obj.comparator.density;
            c(:,1) = C(obj.istre,obj.jstre,:);
            obj.h{1} = plot(rho,c,'+');
        end
        
        function plotVademecumTensor(obj)
            C   = obj.comparator.CtensorVademecum;
            rho = obj.comparator.density;
            c(:,1) = C(obj.istre,obj.jstre,:);
            obj.h{2} = plot(rho,c,'+');
            hold on
        end
        
        function addLegend(obj)
            obj.cleg = ['C_{',num2str(obj.istre),num2str(obj.jstre),'}'];
            leg1 = [obj.cleg,' Vademecum'];
            leg2 = [obj.cleg,' SIMPALL'];
            legend({leg1,leg2})
        end
        
        function printPlot(obj)
            plotName = [obj.cleg,'Comparison'];
            p = plotPrinter(obj.fig,obj.h);
            p.print([obj.path,plotName]);
        end
        
    end
end