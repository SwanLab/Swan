function s = stress_plate(U,p,t,matprop)

    nt=size(t,2); % Number of triangles
    np=size(p,2); % Number of points

    % Corner point indices
    it1=t(1,:);
    it2=t(2,:);
    it3=t(3,:);

    %Triangle coorners coordinates
    %Global coordinates
    x1g=p(1,it1)';
    x2g=p(1,it2)';
    x3g=p(1,it3)';
    y1g=p(2,it1)';
    y2g=p(2,it2)';
    y3g=p(2,it3)';

    %Calcula beta
    beta=zeros(nt,1);
    beta(:) = asin(1.0);
    delx = x2g-x1g;
    dely = y2g-y1g;
    for i=1:nt
        if (delx(i) ~= 0.0)
            beta(i) = atan(dely(i)./delx(i));
        elseif (dely(i) < 0.0)
            beta(i) = -beta(i);
        end
    end

    %Local coordinates
    x1 = cos(beta).*(x1g - x1g)+sin(beta).*(y1g - y1g);
    y1 =-sin(beta).*(x1g - x1g)+cos(beta).*(y1g - y1g);

    x2 = cos(beta).*(x2g - x1g)+sin(beta).*(y2g - y1g);
    y2 =-sin(beta).*(x2g - x1g)+cos(beta).*(y2g - y1g);

    x3 = cos(beta).*(x3g - x1g)+sin(beta).*(y3g - y1g);
    y3 =-sin(beta).*(x3g - x1g)+cos(beta).*(y3g - y1g);

    %Alpha Matrix
    x12 = x1-x2;
    y12 = y1-y2;
    x23 = x2-x3;
    y23 = y2-y3;
    x31 = x3-x1;
    y31 = y3-y1;

    l122 = x12.*x12 + y12.*y12;
    l232 = x23.*x23 + y23.*y23;
    l312 = x31.*x31 + y31.*y31;

    p4 = -6*x23./l232;
    p5 = -6*x3 ./l312;
    p6 = -6*x12./l122;

    t4 = -6*y23./l232;
    t5 = -6*y3 ./l312;

    q4 = 3*x23.*y23./l232;
    q5 = 3*x3.*y3  ./l312;

    r4 = 3*(y23.*y23)./l232;
    r5 = 3*(y31.*y31)./l312;

    TwoArea = x2 .*y3;

    L = 1/3 * [1 1 1 0 0 0 0 0 0;
               0 0 0 1 1 1 0 0 0;
               0 0 0 0 0 0 1 1 1];

    %Material Properties
    h0=matprop.h0;
    nu=matprop.nu0;
    E0=matprop.E0(1);

    Df11=E0*h0*h0*h0/(12*(1-nu*nu));
    Df12=nu*Df11;
    Df22=Df11;
    Df33=Df11*(1-nu)/2;

    Df=[Df11   Df12    0;
        Df12   Df22    0;
        0      0       Df33];

    s=zeros(3,nt);

    for i=1:nt
        %Alfa matrix
        alfa=zeros(9,9);
        alfa(1,1) = y3(i)*p6(i);
        alfa(1,2) = 0.0;
        alfa(1,3) = -4.0*y3(i);
        alfa(1,4) = -y3(i)*p6(i);
        alfa(1,5) = 0.0;
        alfa(1,6) = -2.0*y3(i);
        alfa(1,7) = 0.0;
        alfa(1,8) = 0.0;
        alfa(1,9) = 0.0;
        alfa(2,1) = -y3(i)*p6(i);
        alfa(2,2) = 0.0;
        alfa(2,3) = 2.0*y3(i);
        alfa(2,4) = y3(i)*p6(i);
        alfa(2,5) = 0.0;
        alfa(2,6) = 4.0*y3(i);
        alfa(2,7) = 0.0;
        alfa(2,8) = 0.0;
        alfa(2,9) = 0.0;
        alfa(3,1) = y3(i)*p5(i);
        alfa(3,2) = -y3(i)*q5(i);
        alfa(3,3) = y3(i)*(2.0 - r5(i));
        alfa(3,4) = y3(i)*p4(i);
        alfa(3,5) = y3(i)*q4(i);
        alfa(3,6) = y3(i)*(r4(i) - 2.0);
        alfa(3,7) = -y3(i)*(p4(i) + p5(i));
        alfa(3,8) = y3(i)*(q4(i) - q5(i));
        alfa(3,9) = y3(i)*(r4(i) - r5(i));
        alfa(4,1) = -x2(i)*t5(i);
        alfa(4,2) = x23(i) + x2(i)*r5(i);
        alfa(4,3) = -x2(i)*q5(i);
        alfa(4,4) = 0.0;
        alfa(4,5) = x3(i);
        alfa(4,6) = 0.0;
        alfa(4,7) = x2(i)*t5(i);
        alfa(4,8) = x2(i)*(r5(i) - 1.0);
        alfa(4,9) = -x2(i)*q5(i);
        alfa(5,1) = 0.0;
        alfa(5,2) = x23(i);
        alfa(5,3) = 0.0;
        alfa(5,4) = x2(i)*t4(i);
        alfa(5,5) = x3(i) + x2(i)*r4(i);
        alfa(5,6) = -x2(i)*q4(i);
        alfa(5,7) = -x2(i)*t4(i);
        alfa(5,8) = x2(i)*(r4(i) - 1.0);
        alfa(5,9) = -x2(i)*q4(i);
        alfa(6,1) = x23(i)*t5(i);
        alfa(6,2) = x23(i)*(1.0 - r5(i));
        alfa(6,3) = x23(i)*q5(i);
        alfa(6,4) = -x3(i)*t4(i);
        alfa(6,5) = x3(i)*(1 - r4(i));
        alfa(6,6) = x3(i)*q4(i);
        alfa(6,7) = -x23(i)*t5(i) +x3(i)*t4(i);
        alfa(6,8) = -x23(i)*r5(i) - x3(i)*r4(i) - x2(i);
        alfa(6,9) = x3(i)*q4(i) + x23(i)*q5(i);
        alfa(7,1) = -x3(i)*p6(i) - x2(i)*p5(i);
        alfa(7,2) = x2(i)*q5(i) + y3(i);
        alfa(7,3) = -4.0*x23(i) + x2(i)*r5(i);
        alfa(7,4) = x3(i)*p6(i);
        alfa(7,5) = -y3(i);
        alfa(7,6) = 2.0*x3(i);
        alfa(7,7) = x2(i)*p5(i);
        alfa(7,8) = x2(i)*q5(i);
        alfa(7,9) = (r5(i) - 2.0)*x2(i);
        alfa(8,1) = -x23(i)*p6(i);
        alfa(8,2) = y3(i);
        alfa(8,3) = 2.0*x23(i);
        alfa(8,4) = x23(i)*p6(i) + x2(i)*p4(i);
        alfa(8,5) = -y3(i) + x2(i)*q4(i);
        alfa(8,6) = -4.0*x3(i)+x2(i)*r4(i);
        alfa(8,7) = -x2(i)*p4(i);
        alfa(8,8) = x2(i)*q4(i);
        alfa(8,9) = (r4(i) - 2.0)*x2(i);
        alfa(9,1) = x23(i)*p5(i) + y3(i)*t5(i);
        alfa(9,2) = -x23(i)*q5(i) + (1.0 - r5(i))*y3(i);
        alfa(9,3) = (2.0 - r5(i))*x23(i) + y3(i)*q5(i);
        alfa(9,4) = -x3(i)*p4(i) + y3(i)*t4(i);
        alfa(9,5) = (r4(i) - 1.0)*y3(i) - x3(i)*q4(i);
        alfa(9,6) = (2.0 - r4(i))*x3(i) - y3(i)*q4(i);
        alfa(9,7) = -x23(i)*p5(i) + x3(i)*p4(i) - (t4(i) + t5(i))*y3(i);
        alfa(9,8) = -x23(i)*q5(i) - x3(i)*q4(i) + (r4(i) - r5(i))*y3(i);
        alfa(9,9) = -x23(i)*r5(i) - x3(i)*r4(i) + 4.0*x2(i) + (q5(i) - q4(i))*y3(i);

        %Calculates B
        B=L*alfa;

        T=zeros(9,9);
        T(1,1) =  1;
        T(2,2) =  cos(beta(i));
        T(2,3) =  sin(beta(i));
        T(3,2) = -sin(beta(i));
        T(3,3) =  cos(beta(i));
        T(4,4) =  1;
        T(5,5) =  cos(beta(i));
        T(5,6) =  sin(beta(i));
        T(6,5) = -sin(beta(i));
        T(6,6) =  cos(beta(i));
        T(7,7) =  1;
        T(8,8) =  cos(beta(i));
        T(8,9) =  sin(beta(i));
        T(9,8) = -sin(beta(i));
        T(9,9) =  cos(beta(i));

        %Rotaciona vetor deslocamento [u] = [T][U]
        U1 = zeros(9,1);
        U1(1) = U(it1(i));
        U1(2) = U(it1(i)+np);
        U1(3) = U(it1(i)+2*np);
        U1(4) = U(it2(i));
        U1(5) = U(it2(i)+np);
        U1(6) = U(it2(i)+2*np);
        U1(7) = U(it3(i));
        U1(8) = U(it3(i)+np);
        U1(9) = U(it3(i)+2*np);

        u = T * U1;

        %Calcula deformacao local
        DefL = B * u;

        %Calcula tensao local
        StrL = Df * DefL;

        StrL(1:3) = StrL(1:3)/TwoArea(i);
        %DefL(1:3) = DefL(1:3)/TwoArea(i);

        % Rotaciona Deformacao e tensao
        cc = cos(beta(i))*cos(beta(i));
        ss = sin(beta(i))*sin(beta(i));
        sc = sin(beta(i))*cos(beta(i));

        % Global s = [sxx, sxy, syy]
        s(1,i) = (cc*StrL(1)+ss*StrL(2)- 2.0*sc*StrL(3));
        s(2,i) = (sc*StrL(1)-sc*StrL(2)+(cc-ss)*StrL(3));
        s(3,i) = (ss*StrL(1)+cc*StrL(2)+ 2.0*sc*StrL(3));
    end	
    
end