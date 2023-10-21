classdef plotError < handle 
    
    properties (Access = private)
        residuConjugateB
        residuConjugateRand
        residuConjugateZeros
        residuGaussB
        residuGaussRand
        residuGaussZeros
        residuJacobiB
        residuJacobiRand
        residuJacobiZeros
        
        numIterConjugateB
        numIterConjugateRand
        numIterConjugateZeros
        numIterGaussB
        numIterGaussRand
        numIterGaussZeros
        numIterJacobiB
        numIterJacobiRand
        numIterJacobiZeros
    end
    
    methods (Access = public)
        
        function obj = plotError
            obj.init;
            
            figure
            plot(obj.numIterConjugateB,obj.residuConjugateB.residu,obj.numIterConjugateRand,obj.residuConjugateRand.residu,obj.numIterConjugateZeros,obj.residuConjugateZeros.residu)
            set(findall(gcf,'-property','FontSize'),'FontSize',24)
            set(gca, 'YScale', 'log')
            %set(gca, 'XScale', 'log')
            title('Conjugate Gradient')
            xlabel('Iteration')
            ylabel('Error')
            legend({'B', 'Random', 'Zeros'}, 'Location', 'northeast')
            figure
            plot(obj.numIterGaussB,obj.residuGaussB.residu,obj.numIterGaussRand,obj.residuGaussRand.residu,obj.numIterGaussZeros,obj.residuGaussZeros.residu)
            set(findall(gcf,'-property','FontSize'),'FontSize',24)
            set(gca, 'YScale', 'log')
            %set(gca, 'XScale', 'log')
            title('Gauss-Seidel')
            xlabel('Iteration')
            ylabel('Error')
            legend({'B', 'Random', 'Zeros'}, 'Location', 'northeast')
            %xlim ([0 100])
            figure
            plot(obj.numIterJacobiB,obj.residuJacobiB.residu,obj.numIterJacobiRand,obj.residuJacobiRand.residu,obj.numIterJacobiZeros,obj.residuJacobiZeros.residu)
            set(findall(gcf,'-property','FontSize'),'FontSize',24)
            set(gca, 'YScale', 'log')
            %set(gca, 'XScale', 'log')
            title('Jacobi')
            xlabel('Iteration')
            ylabel('Error')
            legend({'B', 'Random', 'Zeros'}, 'Location', 'northeast')
            %xlim ([0 100])
            
        end
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.residuConjugateB = load('residuConjugateB','residu');
            obj.numIterConjugateB = linspace(1, length(obj.residuConjugateB.residu), length(obj.residuConjugateB.residu));
            obj.residuConjugateRand = load('residuConjugateRand','residu');
            obj.numIterConjugateRand = linspace(1, length(obj.residuConjugateRand.residu), length(obj.residuConjugateRand.residu));
            obj.residuConjugateZeros = load('residuConjugateZeros','residu');
            obj.numIterConjugateZeros = linspace(1, length(obj.residuConjugateZeros.residu), length(obj.residuConjugateZeros.residu));
            obj.residuGaussB = load('residuGaussB','residu');
            obj.numIterGaussB = linspace(1, length(obj.residuGaussB.residu), length(obj.residuGaussB.residu));
            obj.residuGaussRand = load('residuGaussRand','residu');
            obj.numIterGaussRand = linspace(1, length(obj.residuGaussRand.residu), length(obj.residuGaussRand.residu));
            obj.residuGaussZeros = load('residuGaussZeros','residu');
            obj.numIterGaussZeros = linspace(1, length(obj.residuGaussZeros.residu), length(obj.residuGaussZeros.residu));
            obj.residuJacobiB = load('residuJacobiB','residu');
            obj.numIterJacobiB = linspace(1, length(obj.residuJacobiB.residu), length(obj.residuJacobiB.residu));
            obj.residuJacobiRand = load('residuJacobiRand','residu');
            obj.numIterJacobiRand = linspace(1, length(obj.residuJacobiRand.residu), length(obj.residuJacobiRand.residu));
            obj.residuJacobiZeros = load('residuJacobiZeros','residu');
            obj.numIterJacobiZeros = linspace(1, length(obj.residuJacobiZeros.residu), length(obj.residuJacobiZeros.residu));
        end
    end
    
end
