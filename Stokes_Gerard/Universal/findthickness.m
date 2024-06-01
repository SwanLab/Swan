function [t,Q] = findthickness(x_pos,ampl,alt,M,p)

% Trobar t
%Primera posició indica quin node del rectangle és i la segona si és x o y
Q(1,1)=x_pos-ampl/2;
Q(2,1)=x_pos+ampl/2;
Q(3,1)=x_pos+ampl/2;
Q(4,1)=x_pos-ampl/2;

if x_pos<=p
    Q(1,2) = ((M/(p^2))*(2*p*x_pos-x_pos^2) + (alt/2));
    Q(2,2) = ((M/(p^2))*(2*p*x_pos-x_pos^2) + (alt/2));
    Q(3,2) = ((M/(p^2))*(2*p*x_pos-x_pos^2) - (alt/2));
    Q(4,2) = ((M/(p^2))*(2*p*x_pos-x_pos^2) - (alt/2));

else
    Q(1,2) = ((M/(1-p)^2)*((1-2*p)+2*p*x_pos-x_pos^2)) + (alt/2);
    Q(2,2) = ((M/(1-p)^2)*((1-2*p)+2*p*x_pos-x_pos^2)) + (alt/2);
    Q(3,2) = ((M/(1-p)^2)*((1-2*p)+2*p*x_pos-x_pos^2)) - (alt/2);
    Q(4,2) = ((M/(1-p)^2)*((1-2*p)+2*p*x_pos-x_pos^2)) - (alt/2);
end

err = 1*10^(-6);

%Extradós:
for j=1:1:2
    if Q(j,1)<=p
     
        f = @(X) X - Q(j,1) - ((Q(j,2)-((M/(p^2))*(2*p*X-X^2)))/cos(atan((2*M/(p^2))*(p-X))))*sin(atan((2*M/(p^2))*(p-X)));
        A = 0;
        B = p;
        fB = f(B);
        fA = f(A);

        if fA * fB > 0
            error('La funció no canvia de signe en l''interval [A, B].');
        end
        while (abs(A-B)>err)

            C = (A+B)/2;
            fC=f(C);

            if fA*fC > 0
                A = C;
                fA = fC;

            else
                B = C;
                fB = fC;
            end

        end

        theta=atan((2*M/(p^2))*(p-C)); %C serà la x que és l'origen de la circumferència que defineix el punt que considerem del quadrat

        y_t = (Q(j,2)-((M/(p^2))*(2*p*C-C^2)))/cos(theta);
        
        T(j) = y_t/(5*(0.2969*sqrt(C)-0.1260*C-0.3516*C^2+0.2843*C^3-0.1015*C^4));


    elseif Q(j,1)>p
        
        f = @(X) X - Q(j,1) - ((Q(j,2)-((M/(1-p)^2)*((1-2*p)+2*p*X-X^2)))/cos(atan((2*M/((1-p)^2))*(p-X))))*sin(atan((2*M/((1-p)^2))*(p-X)));
        A = p;
        B = 1;
        fB = f(B);
        fA = f(A);

        if fA * fB > 0
            error('La funció no canvia de signe en l''interval [A, B].');
        end
        while (abs(A-B)>err)

            C = (A+B)/2;
            fC=f(C);

            if fA*fC > 0
                A = C;
                fA = fC;

            else
                B = C;
                fB = fC;
            end

        end

        theta=atan(((2*M/((1-p)^2))*(p-C))); %C serà la x que és l'origen de la circumferència que defineix el punt que considerem del quadrat

        y_t = (Q(j,2)-((M/(1-p)^2)*((1-2*p)+2*p*C-C^2)))/cos(theta);
        
        T(j) = y_t/(5*(0.2969*sqrt(C)-0.1260*C-0.3516*C^2+0.2843*C^3-0.1015*C^4));
    end

end


%Intradós:

for j=3:1:4
    if Q(j,1)<=p
             
        f = @(X) X - Q(j,1) - ((Q(j,2)-((M/(p^2))*(2*p*X-X^2)))/cos(atan((2*M/(p^2))*(p-X))))*sin(atan((2*M/(p^2))*(p-X)));
        A = 0;
        B = p;
        fB = f(B);
        fA = f(A);

        if fA * fB > 0
            error('La funció no canvia de signe en l''interval [A, B].');
        end
        while (abs(A-B)>err)

            C = (A+B)/2;
            fC=f(C);

            if fA*fC > 0
                A = C;
                fA = fC;

            else
                B = C;
                fB = fC;
            end

        end

        theta=atan((2*M/(p^2))*(p-C)); %C serà la x que és l'origen de la circumferència que defineix el punt que considerem del quadrat

        y_t = (-Q(j,2)+((M/(p^2))*(2*p*C-C^2)))/cos(theta);
        
        T(j) = y_t/(5*(0.2969*sqrt(C)-0.1260*C-0.3516*C^2+0.2843*C^3-0.1015*C^4));


    elseif Q(j,1)>p
        f = @(X) X - Q(j,1) - ((Q(j,2)-((M/(1-p)^2)*((1-2*p)+2*p*X-X^2)))/cos(atan((2*M/((1-p)^2))*(p-X))))*sin(atan((2*M/((1-p)^2))*(p-X)));
        A = p;
        B = 1;
        fB = f(B);
        fA = f(A);

        if fA * fB > 0
            error('La funció no canvia de signe en l''interval [A, B].');
        end
        while (abs(A-B)>err)

            C = (A+B)/2;
            fC=f(C);

            if fA*fC > 0
                A = C;
                fA = fC;

            else
                B = C;
                fB = fC;
            end

        end

        theta=atan(((2*M/((1-p)^2))*(p-C))); %C serà la x que és l'origen de la circumferència que defineix el punt que considerem del quadrat

        y_t = (-Q(j,2)+((M/(1-p)^2)*((1-2*p)+2*p*C-C^2)))/cos(theta);
        
        T(j) = y_t/(5*(0.2969*sqrt(C)-0.1260*C-0.3516*C^2+0.2843*C^3-0.1015*C^4));


    end

end


[t, maxIndex_t] = max(T);


% hold on
% scatter(Q(:,1),Q(:,2))

end












































