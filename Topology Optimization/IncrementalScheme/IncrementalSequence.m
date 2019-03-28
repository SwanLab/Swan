classdef IncrementalSequence < handle
    
    properties (GetAccess = public, SetAccess = private)
        value
    end
    
    properties (Access = public)
        alpha
        initialValue
        finalValue     
    end
    
    methods (Access = public)
        
        function obj = IncrementalSequence(x0,x1,nsteps,type,initialValue,finalValue,factor)
            obj.init(initialValue,finalValue);
           obj.generateAlphaSequence(x0,x1,nsteps,type); 
        end
        
        function update(obj,i)
            obj.value = (1-obj.alpha(i))*obj.initialValue + obj.alpha(i)*obj.finalValue;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,a0,a1)
            obj.initialValue = a0;
            obj.finalValue = a1;
        end
        
        function generateAlphaSequence(obj,x0,x1,nsteps,type,factor)
            switch type
                case 'linear'
                    x = linspace(x0,x1,nsteps);
                case 'epsilon_sequence'
                    frac = 2;
                    kmax = ceil(log10(x0/x1)/log10(frac));
                    x = obj.epsilon0./frac.^(1:kmax);
                case 'logarithmic'
                    x = logspace(x0,x1,nsteps);
                case 'custom'
                    if nsteps < 2
                        x = x1;
                    else
                        isteps = 0:nsteps-1;
                        x = 1-(1-isteps/(nsteps-1)).^(factor);
                        x = (x1-x0)*x + x0;
                    end
                case 'free'
                    x = zeros(1,nsteps);
                    x(end) = 1;
                otherwise
                    error('Incremental sequence type not detected.')
            end
            obj.alpha = x;
        end
        
    end
    
end

