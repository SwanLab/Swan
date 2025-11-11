function fittingWithSplines()
    close all
    % % C12data = squeeze(matType.mat(1,1,2,2,:));
    % % C33data = squeeze(matType.mat(1,2,1,2,:));
    % % Cdata = [C11data';C12data';C33data'];
    % % 
    % % xq = 0:1e-3:1;
    % % fun = griddedInterpolant([phiData 1],[C11data' 0],'spline');
    % % vq = fun(xq);
    % % plot(xq,vq);

    matType = load('/home/gerard/Documents/GitHub/Swan/src/Problems/Damage/Models/PhaseField/PFVademecum/Degradation/SquareArea.mat');

    x = matType.phi;
    y = squeeze(matType.mat(1,1,1,1,:));


    Alinear = [1 x(1) x(1)^2 x(1)^3;
               1 x(2) x(2)^2 x(2)^3];
    Blinear = [y(1); y(2)];
    for i=2:length(x)-1
        values = [1 x(i)   x(i)^2   x(i)^3;
                  1 x(i+1) x(i+1)^2 x(i+1)^3];
        Alinear = [     Alinear             zeros(size(Alinear,1),4);
                    zeros(2,size(Alinear,2))      values        ];
        Blinear = [Blinear; y(i); y(i+1)];
    end

    Adiff = [0  1  2*x(1)  3*x(1)^2;
             0 -1 -2*x(2) -3*x(2)^2];
    for i=2:length(x)-1
        values = [0  1  2*x(i)      3*x(i)^2;
                  0 -1 -2*x(i+1) -3*x(i+1)^2];
        Adiff = [[Adiff;zeros(1,size(Adiff,2))] [zeros(size(Adiff,1)-1,4);values]];
    end
    k = -5;
    Bdiff = [k;zeros(size(Adiff,1)-1,1)];

    A2diff = [0 0 2 6*x(2) 0  0 -2 -6*x(2)];
    for i=3:length(x)-1
        values = [0 0 2 6*x(i) 0  0 -2 -6*x(i)];
        mat = zeros(size(A2diff,1)+1,size(A2diff,2)+4);
        mat(1:end-1,1:end-4) = A2diff;
        mat(end,end-7:end) = values;
        A2diff = mat;
    end
    B2diff = zeros(size(A2diff,1),1);

    A = [Alinear;Adiff;A2diff];
    B = [Blinear;Bdiff;B2diff];

    coeff = A\B;
    
    syms xFun
    spline = cell(1,length(x)-1);
    for i=1:length(x)-1
        init = 4*(i-1)+1;
        spline{i} = coeff(init) + coeff(init+1)*xFun + coeff(init+2)*xFun^2 + coeff(init+3)*xFun^3;
    end

    figure(1)
    hold on
    for i=1:length(x)-1
        fun = spline{i};
        fplot(fun,[x(i) x(i+1)]);
    end

    figure(2)
    hold on
    for i=1:length(x)-1
        fun = diff(spline{i});
        fplot(fun,[x(i) x(i+1)]);
    end

    figure(3)
    hold on
    for i=1:length(x)-1
        fun = diff(diff(spline{i}));
        fplot(fun,[x(i) x(i+1)]);
    end

end