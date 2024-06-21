function f = SourceTerm(x, type, t) 
    % Only time dependant source term
    switch type
        case 'none'
            f = 0;
        case 'unitary'
            f = 1;
        case 'expDecay'
            A = 1;
            alpha = 0.5;
            f = A*exp(-alpha*t);
        case 'step'
            A = 1;
            t1 = 0;
            t2 = 2;
            if t >= t1 && t <= t2
                f = A;
            else
                f = 0;
            end
    end
end
