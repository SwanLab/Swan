classdef OrderTextInterpreter 
   
    methods (Static,Access = public)
        function ord = compute(txt)
            switch (txt)
                case 'LINEAR'
                    ord = 1;
                case 'QUADRATIC'
                    ord = 2;
                case 'CUBIC'
                    ord = 3;
            end
        end
    end
    
end