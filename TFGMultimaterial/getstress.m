% stress computation for PLATES
function s = getstress(U,mesh,matprop,pdecoef)

p = mesh.p; t = mesh.t;
sf = size(pdecoef.f,1);
sc = size(pdecoef.c,1);
nt=size(t,2); % Number of triangles
np=size(p,2); % Number of points
    
if sf==3 % 
    if sc==2 % stress computation for Kirchhoff plates
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

        auxdelx = (delx ~= 0.0);
        auxdely = (dely < 0.0);
        beta(auxdelx) = atan(dely(auxdelx)./delx(auxdelx));
        beta(~auxdelx & auxdely) = -beta(~auxdelx & auxdely);
        Cos = cos(beta);
        Sen = sin(beta);

        %Local coordinates
        x1 = Cos.*(x1g - x1g) + Sen.*(y1g - y1g);
        y1 =-Sen.*(x1g - x1g) + Cos.*(y1g - y1g);

        x2 = Cos.*(x2g - x1g) + Sen.*(y2g - y1g);
        y2 =-Sen.*(x2g - x1g) + Cos.*(y2g - y1g);

        x3 = Cos.*(x3g - x1g) + Sen.*(y3g - y1g);
        y3 =-Sen.*(x3g - x1g) + Cos.*(y3g - y1g);

        x1=x1';
        y1=y1';
        x2=x2';
        y2=y2';
        x3=x3';
        y3=y3';
        Cos = Cos';
        Sen = Sen';

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

        %Material Properties

        h0 = matprop.h0; 
%         la = matprop.la; mu = matprop.mu;     
%         nu = la/(2*mu+la); E0 = 2*mu*(1+nu);
        nu = matprop.nu;
        E0 = matprop.E(1);
        
        %Material Properties    
        E1 = E0*h0*h0*h0/(12*(1-nu*nu));
        E2 = nu.*E1;
        E3 = E1;
        E4 = E1.*(1.-nu)./2;


        %Elements of stiffness matrix
        alfa11 = y3.*p6;
        %alfa12 = 0.0;
        alfa13 = -4.0.*y3;
        alfa14 = -alfa11;
        %alfa(1,5) = 0.0;
        alfa16 = -2.0.*y3;
        %alfa17 = 0.0;
        %alfa18 = 0.0;
        %alfa19 = 0.0;
        alfa21 = alfa14;
        %alfa22 = 0.0;
        alfa23 = -alfa16;
        alfa24 = alfa11;
        %alfa25 = 0.0;
        alfa26 = -alfa13;
        %alfa27 = 0.0;
        %alfa28 = 0.0;
        %alfa29 = 0.0;
        alfa31 = y3.*p5;
        alfa32 = -y3.*q5;
        alfa33 = y3.*(2.0 - r5);
        alfa34 = y3.*p4;
        alfa35 = y3.*q4;
        alfa36 = y3.*(r4 - 2.0);
        alfa37 = -y3.*(p4 + p5);
        alfa38 = y3.*(q4 - q5);
        alfa39 = y3.*(r4 - r5);
        alfa41 = -x2.*t5;
        alfa42 = x23 + x2.*r5;
        alfa43 = -x2.*q5;
        %alfa44 = 0.0;
        alfa45 = x3;
        %alfa46 = 0.0;
        alfa47 = -alfa41;
        alfa48 = x2.*(r5 - 1.0);
        alfa49 = alfa43;
        %alfa51 = 0.0;
        alfa52 = x23;
        %alfa53 = 0.0;
        alfa54 = x2.*t4;
        alfa55 = x3 + x2.*r4;
        alfa56 = -x2.*q4;
        alfa57 = -alfa54;
        alfa58 = x2.*(r4 - 1.0);
        alfa59 = alfa56;
        alfa61 = x23.*t5;
        alfa62 = x23.*(1.0 - r5);
        alfa63 = x23.*q5;
        alfa64 = -x3.*t4;
        alfa65 = x3.*(1.0 - r4);
        alfa66 = x3.*q4;
        alfa67 = -alfa61 -alfa64;
        alfa68 = -x23.*r5 - x3.*r4 - x2;
        alfa69 = alfa66 + alfa63;
        alfa71 = -x3.*p6 - x2.*p5;
        alfa72 = -alfa43 + y3;
        alfa73 = -4.0.*x23 + x2.*r5;
        alfa74 = x3.*p6;
        alfa75 = -y3;
        alfa76 = 2.0.*x3;
        alfa77 = x2.*p5;
        alfa78 = -alfa43;
        alfa79 = (r5 - 2.0).*x2;
        alfa81 = -x23.*p6;
        alfa82 = y3;
        alfa83 = 2.0.*x23;
        alfa84 = -alfa81 + x2.*p4;
        alfa85 = -y3 - alfa56;
        alfa86 = -4.0.*x3 + x2.*r4;
        alfa87 = -x2.*p4;
        alfa88 = -alfa56;
        alfa89 = (r4 - 2.0).*x2;
        alfa91 = x23.*p5 + y3.*t5;
        alfa92 = -alfa63 + (1.0 - r5).*y3;
        alfa93 = (2.0 - r5).*x23 - alfa32;
        alfa94 = -x3.*p4 + y3.*t4;
        alfa95 = (r4 - 1.0).*y3 - alfa66;
        alfa96 = (2.0 - r4).*x3 - alfa35;
        alfa97 = -x23.*p5 + x3.*p4 - (t4 + t5).*y3;
        alfa98 = -alfa63 - alfa66 + (r4 - r5).*y3;
        alfa99 = -x23.*r5 - x3.*r4 + 4.0*x2 + (q5 - q4).*y3;

        % vetor deslocamento
        U1 = U(it1)';
        U2 = U(it1+np)';
        U3 = U(it1+2*np)';
        U4 = U(it2)';
        U5 = U(it2+np)';
        U6 = U(it2+2*np)';
        U7 = U(it3)';
        U8 = U(it3+np)';
        U9 = U(it3+2*np)';

        StrL(1,:) = (1/3./TwoArea).*((((alfa11 + alfa21 + alfa31).*U1 + (alfa14 + alfa24 + alfa34).*U4 + alfa37.*U7) ...
                  + (alfa32.*U2 + (alfa13 + alfa23 + alfa33).*U3 + alfa35.*U5 + (alfa16 + alfa26 + alfa36).*U6 + alfa38.*U8 + alfa39.*U9).*Cos - ((alfa13 + alfa23 + alfa33).*U2 - alfa32.*U3 + (alfa16 + alfa26 + alfa36).*U5 - alfa35.*U6 + alfa39.*U8 - alfa38.*U9).*Sen).*E1 ...
                  + (((alfa41 + alfa61).*U1 + (alfa54 + alfa64).*U4 + (alfa47 + alfa57 + alfa67).*U7) ...
                  + ((alfa42 + alfa52 + alfa62).*U2 + (alfa43 + alfa63).*U3 + (alfa45 + alfa55 + alfa65).*U5 + (alfa56 + alfa66).*U6 + (alfa48 + alfa58 + alfa68).*U8 + (alfa49 + alfa59 + alfa69).*U9).*Cos - ((alfa43 + alfa63).*U2 - (alfa42 + alfa52 + alfa62).*U3 + (alfa56 + alfa66).*U5 - (alfa45 + alfa55 + alfa65).*U6 + (alfa49 + alfa59 + alfa69).*U8 - (alfa48 + alfa58 + alfa68).*U9).*Sen).*E2);

        StrL(2,:) = (1/3./TwoArea).*((((alfa11 + alfa21 + alfa31).*U1 + (alfa14 + alfa24 + alfa34).*U4 + alfa37.*U7) ...
                  + (alfa32.*U2 + (alfa13 + alfa23 + alfa33).*U3 + alfa35.*U5 + (alfa16 + alfa26 + alfa36).*U6 + alfa38.*U8 + alfa39.*U9).*Cos - ((alfa13 + alfa23 + alfa33).*U2 - alfa32.*U3 + (alfa16 + alfa26 +alfa36).*U5 - alfa35.*U6 + alfa39.*U8 - alfa38.*U9).*Sen).*E2 ...
                  + (((alfa41 + alfa61).*U1 + (alfa54 + alfa64).*U4 + (alfa47 + alfa57 + alfa67).*U7) ...
                  + ((alfa42 + alfa52 + alfa62).*U2 + (alfa43 + alfa63).*U3 + (alfa45 + alfa55 + alfa65).*U5 + (alfa56 + alfa66).*U6 + (alfa48 + alfa58 + alfa68).*U8 + (alfa49 + alfa59 + alfa69).*U9).*Cos - ((alfa43 + alfa63).*U2 - (alfa42 + alfa52 + alfa62).*U3 + (alfa56 + alfa66).*U5 - (alfa45 + alfa55 + alfa65).*U6 + (alfa49 + alfa59 + alfa69).*U8 - (alfa48 + alfa58 + alfa68).*U9).*Sen).*E3);

        StrL(3,:) = (1/3./TwoArea).*((((alfa71 + alfa81 + alfa91).*U1 + (alfa74 + alfa84 + alfa94).*U4 + (alfa77 + alfa87 + alfa97).*U7) ...
                  + ((alfa72 + alfa82 + alfa92).*U2 + (alfa73 + alfa83 + alfa93).*U3 + (alfa75 + alfa85 + alfa95).*U5 + (alfa76 + alfa86 + alfa96).*U6 + (alfa78 + alfa88 + alfa98).*U8 + (alfa79 + alfa89 + alfa99).*U9).*Cos - ((alfa73 + alfa83 + alfa93).*U2 - (alfa72 + alfa82 + alfa92).*U3 + (alfa76 + alfa86 + alfa96).*U5 - (alfa75 + alfa85 + alfa95).*U6 + (alfa79 + alfa89 + alfa99).*U8 - (alfa78 + alfa88 + alfa98).*U9).*Sen).*E4);

        % Rotaciona Deformacao e tensao
        cc = Cos.*Cos;
        ss = Sen.*Sen;
        sc = Sen.*Cos;

        % Global s = [sxx, sxy, syy]
        s=zeros(3,nt);
        s(1,:) = (cc.*StrL(1,:) + ss.*StrL(2,:) - 2.0.*sc.*StrL(3,:));
        s(2,:) = (sc.*StrL(1,:) - sc.*StrL(2,:) + (cc - ss).*StrL(3,:));
        s(3,:) = (ss.*StrL(1,:) + cc.*StrL(2,:) + 2.0.*sc.*StrL(3,:));
        
    elseif sc==3 % stress computation for Reissner-Mindlin plates
        %============================================================
        % Finite Element implementation based on the work:
        % A NEW DISCRETE KIRCHHOFF-MINDLIN ELEMENT BASED ON MINDLIN-REISSNER
        % PLATE THEORY AND ASSUMED SHEAR STRAIN FIELDS-PART I: AN EXTENDED
        % DKT ELEMENT FOR THICK-PLATE BENDING ANALYSIS
        % by Katili I.
        % in International JOurnal for Numerical Methods in Engineering 36(1993):
        % 1859-1883
        %============================================================
        
        % Corner point indices
        it1=t(1,:);
        it2=t(2,:);
        it3=t(3,:);

        % Triangle geometries:
        ar = mesh.area;

        % Find midpoints of triangles
        xc=(p(1,it1)+p(1,it2)+p(1,it3))/3;
        yc=(p(2,it1)+p(2,it2)+p(2,it3))/3;

        %Triangle coorners coordinates
        %Global coordinates
        x1g=p(1,it1);
        x2g=p(1,it2);
        x3g=p(1,it3);
        y1g=p(2,it1);
        y2g=p(2,it2);
        y3g=p(2,it3);

        % Auxiliary parameters
        x13 = x1g - x3g;
        x21 = x2g - x1g;
        x32 = x3g - x2g;

        y13 = y1g - y3g;
        y21 = y2g - y1g;
        y32 = y3g - y2g;

        L4 = (x21.^2 + y21.^2).^(0.5);
        L5 = (x32.^2 + y32.^2).^(0.5);
        L6 = (x13.^2 + y13.^2).^(0.5);

        C4 = x21./L4;
        C5 = x32./L5;
        C6 = x13./L6;

        S4 = y21./L4;
        S5 = y32./L5;
        S6 = y13./L6;

        L24 = L4.^2;
        L25 = L5.^2;
        L26 = L6.^2;

        %Material Properties
        h0=matprop.h0;
        nu=matprop.nu0;
        E0=matprop.E0(1);
        k = 5/6;

        %Material Properties    
        Ebb = E0*h0*h0*h0/(12*(1-nu*nu));
        Ess = k*E0*h0/(2*(1+nu));
        Z = Ebb./Ess.*6;

        % Shape Function on mid-nodes

        a1=x2g.*y3g-x3g.*y2g;
        a2=x3g.*y1g-x1g.*y3g;
        a3=x1g.*y2g-x2g.*y1g;

        b1=y2g-y3g;
        b2=y3g-y1g;
        b3=y1g-y2g;

        c1=x3g-x2g;
        c2=x1g-x3g;
        c3=x2g-x1g;

        s1 = 1./(2.*ar).*(a1+b1.*xc+c1.*yc);
        s2 = 1./(2.*ar).*(a2+b2.*xc+c2.*yc);
        s3 = 1./(2.*ar).*(a3+b3.*xc+c3.*yc);

        %Elements of Bending matrix

        Mb11 = Ebb.*nu.*((3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(4.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*((3.*C4.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C6.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(4.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb12 = Ebb.*nu.*((-3.*S4.*x21.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S6.*x13.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*(-y32./(2.*ar) - (3.*C4.*x21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C6.*x13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb13 = Ebb.*nu.*(x32./(2.*ar) - (3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21).*y21)./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21).*y13)./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*((-3.*C4.*y21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C6.*y13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb14 = Ebb.*nu.*((-3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) + (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21))./(4.*ar.*L5.*(1 + (2.*Z)./L25))) + Ebb.*((-3.*C4.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) + (3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21))./(4.*ar.*L5.*(1 + (2.*Z)./L25)));
        Mb15 = Ebb.*nu.*((-3.*S4.*x21.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*x32)./(8.*ar.*L5.*(1 + (2.*Z)./L25))) + Ebb.*(-y13./(2.*ar) - (3.*C4.*x21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C5.*x32.*(-4.*s2.*y13 - 4.*s1.*y21))./(8.*ar.*L5.*(1 + (2.*Z)./L25)));
        Mb16 = Ebb.*nu.*(x13./(2.*ar) - (3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21).*y21)./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25))) + Ebb.*((-3.*C4.*y21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)));
        Mb17 = Ebb.*nu.*((-3.*S5.*(4.*s2.*x13 + 4.*s1.*x21))./(4.*ar.*L5.*(1 + (2.*Z)./L25)) + (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(4.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*((-3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21))./(4.*ar.*L5.*(1 + (2.*Z)./L25)) + (3.*C6.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(4.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb18 = Ebb.*nu.*((-3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*x32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*S6.*x13.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*(-y21./(2.*ar) - (3.*C5.*x32.*(-4.*s2.*y13 - 4.*s1.*y21))./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*C6.*x13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb19 = Ebb.*nu.*(x21./(2.*ar) - (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21).*y13)./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*((-3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*C6.*y13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));

        Mb21 = Ebb.*((3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(4.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*nu.*((3.*C4.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C6.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(4.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb22 = Ebb.*((-3.*S4.*x21.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S6.*x13.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*nu.*(-y32./(2.*ar) - (3.*C4.*x21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C6.*x13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb23 = Ebb.*(x32./(2.*ar) - (3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21).*y21)./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21).*y13)./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*nu.*((-3.*C4.*y21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C6.*y13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb24 = Ebb.*((-3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) + (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21))./(4.*ar.*L5.*(1 + (2.*Z)./L25))) + Ebb.*nu.*((-3.*C4.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) + (3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21))./(4.*ar.*L5.*(1 + (2.*Z)./L25)));
        Mb25 = Ebb.*((-3.*S4.*x21.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*x32)./(8.*ar.*L5.*(1 + (2.*Z)./L25))) + Ebb.*nu.*(-y13./(2.*ar) - (3.*C4.*x21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C5.*x32.*(-4.*s2.*y13 - 4.*s1.*y21))./(8.*ar.*L5.*(1 + (2.*Z)./L25)));
        Mb26 = Ebb.*(x13./(2.*ar) - (3.*S4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21).*y21)./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25))) + Ebb.*nu.*((-3.*C4.*y21.*(-4.*(-s1 + s3).*y13 + 4.*s1.*y21))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)));
        Mb27 = Ebb.*((-3.*S5.*(4.*s2.*x13 + 4.*s1.*x21))./(4.*ar.*L5.*(1 + (2.*Z)./L25)) + (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(4.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*nu.*((-3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21))./(4.*ar.*L5.*(1 + (2.*Z)./L25)) + (3.*C6.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(4.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb28 = Ebb.*((-3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*x32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*S6.*x13.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21))./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*nu.*(-y21./(2.*ar) - (3.*C5.*x32.*(-4.*s2.*y13 - 4.*s1.*y21))./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*C6.*x13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));
        Mb29 = Ebb.*(x21./(2.*ar) - (3.*S5.*(4.*s2.*x13 + 4.*s1.*x21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*S6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21).*y13)./(8.*ar.*L6.*(1 + (2.*Z)./L26))) + Ebb.*nu.*((-3.*C5.*(-4.*s2.*y13 - 4.*s1.*y21).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*C6.*y13.*(4.*s2.*y13 - 4.*(-s2 + s3).*y21))./(8.*ar.*L6.*(1 + (2.*Z)./L26)));

        Mb31 = (Ebb.*(1 - nu).*((3.*(C4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21) - S4.*(4.*(-s1 + s3).*y13 - 4.*s1.*y21)))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*(C6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21) - S6.*(-4.*s2.*y13 + 4.*(-s2 + s3).*y21)))./(4.*ar.*L6.*(1 + (2.*Z)./L26))))./2;
        Mb32 = (Ebb.*(1 - nu).*(x32./(2.*ar) - (3.*x21.*(C4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21) - S4.*(4.*(-s1 + s3).*y13 - 4.*s1.*y21)))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*x13.*(C6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21) - S6.*(-4.*s2.*y13 + 4.*(-s2 + s3).*y21)))./(8.*ar.*L6.*(1 + (2.*Z)./L26))))./2;
        Mb33 = (Ebb.*(1 - nu).*(-y32./(2.*ar) - (3.*y21.*(C4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21) - S4.*(4.*(-s1 + s3).*y13 - 4.*s1.*y21)))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*y13.*(C6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21) - S6.*(-4.*s2.*y13 + 4.*(-s2 + s3).*y21)))./(8.*ar.*L6.*(1 + (2.*Z)./L26))))./2;
        Mb34 = (Ebb.*(1 - nu).*((-3.*(C4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21) - S4.*(4.*(-s1 + s3).*y13 - 4.*s1.*y21)))./(4.*ar.*L4.*(1 + (2.*Z)./L24)) + (3.*(C5.*(4.*s2.*x13 + 4.*s1.*x21) - S5.*(4.*s2.*y13 + 4.*s1.*y21)))./(4.*ar.*L5.*(1 + (2.*Z)./L25))))./2;
        Mb35 = (Ebb.*(1 - nu).*(x13./(2.*ar) - (3.*x21.*(C4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21) - S4.*(4.*(-s1 + s3).*y13 - 4.*s1.*y21)))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*x32.*(C5.*(4.*s2.*x13 + 4.*s1.*x21) - S5.*(4.*s2.*y13 + 4.*s1.*y21)))./(8.*ar.*L5.*(1 + (2.*Z)./L25))))./2;
        Mb36 = (Ebb.*(1 - nu).*(-y13./(2.*ar) - (3.*y21.*(C4.*(4.*(-s1 + s3).*x13 - 4.*s1.*x21) - S4.*(4.*(-s1 + s3).*y13 - 4.*s1.*y21)))./(8.*ar.*L4.*(1 + (2.*Z)./L24)) - (3.*(C5.*(4.*s2.*x13 + 4.*s1.*x21) - S5.*(4.*s2.*y13 + 4.*s1.*y21)).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25))))./2;
        Mb37 = (Ebb.*(1 - nu).*((-3.*(C5.*(4.*s2.*x13 + 4.*s1.*x21) - S5.*(4.*s2.*y13 + 4.*s1.*y21)))./(4.*ar.*L5.*(1 + (2.*Z)./L25)) + (3.*(C6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21) - S6.*(-4.*s2.*y13 + 4.*(-s2 + s3).*y21)))./(4.*ar.*L6.*(1 + (2.*Z)./L26))))./2;
        Mb38 = (Ebb.*(1 - nu).*(x21./(2.*ar) - (3.*x32.*(C5.*(4.*s2.*x13 + 4.*s1.*x21) - S5.*(4.*s2.*y13 + 4.*s1.*y21)))./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*x13.*(C6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21) - S6.*(-4.*s2.*y13 + 4.*(-s2 + s3).*y21)))./(8.*ar.*L6.*(1 + (2.*Z)./L26))))./2;
        Mb39 = (Ebb.*(1 - nu).*(-y21./(2.*ar) - (3.*(C5.*(4.*s2.*x13 + 4.*s1.*x21) - S5.*(4.*s2.*y13 + 4.*s1.*y21)).*y32)./(8.*ar.*L5.*(1 + (2.*Z)./L25)) - (3.*y13.*(C6.*(-4.*s2.*x13 + 4.*(-s2 + s3).*x21) - S6.*(-4.*s2.*y13 + 4.*(-s2 + s3).*y21)))./(8.*ar.*L6.*(1 + (2.*Z)./L26))))./2;

        %Elements of Shear matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Original implementation from Katili paper
        
%         [Wx, Wy] = pdegrad(p,t,U(1:np));
%         Gx = Wx + pdeintrp(p,t,U(np+1:2*np));
%         Gy = Wy + pdeintrp(p,t,U(2*np+1:3*np));

%         % vetor deslocamento
%         U1 = U(it1)';
%         U2 = U(it1+np)';
%         U3 = U(it1+2*np)';
%         U4 = U(it2)';
%         U5 = U(it2+np)';
%         U6 = U(it2+2*np)';
%         U7 = U(it3)';
%         U8 = U(it3+np)';
%         U9 = U(it3+2*np)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [Wx, Wy] = pdegrad(p,t,U(1:np));
        Gx = Wy + pdeintrp(p,t,U(np+1:2*np));
        Gy = Wx + pdeintrp(p,t,U(2*np+1:3*np));
        
        % vetor deslocamento
        U1 = U(it1)';
        U3 = U(it1+np)';
        U2 = U(it1+2*np)';
        U4 = U(it2)';
        U6 = U(it2+np)';
        U5 = U(it2+2*np)';
        U7 = U(it3)';
        U9 = U(it3+np)';
        U8 = U(it3+2*np)';

        
        
        % Global s = [Mxx, Myy, Mxy, Sx, Sy]
        s=zeros(5,nt);
        s(1,:) = Mb11.*U1+Mb12.*U2+Mb13.*U3+Mb14.*U4+Mb15.*U5+Mb16.*U6+Mb17.*U7+Mb18.*U8+Mb19.*U9;
        s(2,:) = Mb21.*U1+Mb22.*U2+Mb23.*U3+Mb24.*U4+Mb25.*U5+Mb26.*U6+Mb27.*U7+Mb28.*U8+Mb29.*U9;
        s(3,:) = Mb31.*U1+Mb32.*U2+Mb33.*U3+Mb34.*U4+Mb35.*U5+Mb36.*U6+Mb37.*U7+Mb38.*U8+Mb39.*U9;
        s(4,:) = Ess.*Gx;
        s(5,:) = Ess.*Gy;
                
    end
    
elseif sf==5 % stress computation for Reissner-Mindlin plates
    
    %============================================================
    % Finite Element implementation based on the work:
    % A new and simple locking-free triangular thick plate element
    % using independent shear degrees of freedom
    % by Zhuang, X.Y.; Huang, R.Q.; Zhu, H.H.; Askes, H. & Mathisen, K.
    % in Finite Elements in Analysis and Design 75(2013): 1-7
    %============================================================
    
    % Corner point indices
    it1=t(1,:);
    it2=t(2,:);
    it3=t(3,:);
    
    % Triangle geometries:
    %[ar,~,~,~,~,~,~]=pdetrg(p,t);
    ar = mesh.area;

    % Find midpoints of triangles
    xx=(p(1,it1)+p(1,it2)+p(1,it3))/3;
    yy=(p(2,it1)+p(2,it2)+p(2,it3))/3;
    
    %Triangle coorners coordinates
    %Global coordinates
    x1=p(1,it1);
    y1=p(2,it1);
    x2=p(1,it2);
    y2=p(2,it2);
    x3=p(1,it3);
    y3=p(2,it3);
    
 
    %Material Properties
    h0=matprop.h0;
    nu=matprop.nu0;
    E0=matprop.E0;
    k = 5/6;
  
    Db = E0*h0*h0*h0/(12*(1-nu*nu));
    Ds = k*E0*h0/(2*(1+nu));
    
    %Elements of Bending matrix
    Mb11=-1./2.*Db.*nu.*(x2 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar - 1./2.*Db.*(y2 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar;	
    Mb12=1./8.*Db.*nu.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(2.*x1 - x2 - x3).*(x1 - x2).*(x2 - x3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(2.*x1 - x2 - x3).*(x1 - x3).*(x2 - x3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(x2 - x3).^2./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*(x2 - x3).^2./ar.^3) + 1./8.*Db.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).*((2.*x1 - x2 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).*((2.*x1 - x2 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3).*((x1 - x2).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).*((x1 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 - x2).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y2 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y2 - y3)./ar.^2 + 4.*(y2 - y3).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy)./ar - 1)./ar + 4.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y2 - y3)./ar.^2);	
    Mb13=1./8.*Db.*nu.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x2 - x3).*(2.*y1 - y2 - y3)./ar - ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x2 - x3).*(2.*y1 - y2 - y3)./ar - ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x2 - x3).*(y1 - y2)./ar - ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x2 - x3).*(y1 - y2)./ar - ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x2 - x3).*(y1 - y3)./ar - ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x2 - x3).*(y1 - y3)./ar - ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - 4.*(x2 - x3).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy)./ar - 1)./ar - 4.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x2 - x3)./ar.^2) + 1./8.*Db.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(2.*y1 - y2 - y3).*(y1 - y2).*(y2 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(2.*y1 - y2 - y3).*(y1 - y3).*(y2 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).*(y2 - y3).^2./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).*(y2 - y3).^2./ar.^3);	
    Mb14=-1./2.*Db.*(y2 - y3)./ar;	
    Mb15=1./2.*Db.*nu.*(x2 - x3)./ar;	
    Mb16=1./2.*Db.*nu.*(x1 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar + 1./2.*Db.*(y1 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar;	
    Mb17=-1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(x1 - 2.*x2 + x3).*(x1 - x3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(x1 - x3).^2./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - 2.*x2 + x3).*(x1 - x3).*(x2 - x3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).^2.*(x2 - x3)./ar.^3) - 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - 2.*x2 + x3).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y1 - y2)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x2 - x3).*(y1 - y3)./ar - ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x2).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y1 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x2 - x3).*(y1 - y3)./ar - ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar).*(y1 - y3)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 - x2).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - 2.*x2 + x3).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y2 - y3)./ar.^2 - 4.*(y1 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy)./ar + 1)./ar - 4.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y3)./ar.^2);	
    Mb18=1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x1 - x3).*(y1 - y2)./ar + ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x1 - x3).*(y1 - y2)./ar + ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x1 - x3).*(y1 - 2.*y2 + y3)./ar + ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x1 - x3).*(y1 - 2.*y2 + y3)./ar + ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x1 - x3).*(y2 - y3)./ar + ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x1 - x3).*(y2 - y3)./ar + ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 - 4.*(x1 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy)./ar + 1)./ar - 4.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x3)./ar.^2) - 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).*(y1 - 2.*y2 + y3).*(y1 - y3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).*(y1 - y3).^2./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - 2.*y2 + y3).*(y1 - y3).*(y2 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).^2.*(y2 - y3)./ar.^3);	
    Mb19=1./2.*Db.*(y1 - y3)./ar;	
    Mb110=-1./2.*Db.*nu.*(x1 - x3)./ar;	
    Mb111=-1./2.*Db.*nu.*(x1 - x2).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar - 1./2.*Db.*(y1 - y2).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar;	
    Mb112=-1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 + x2 - 2.*x3).*(x1 - x2).*(x1 - x3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).^2.*(x1 - x3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 + x2 - 2.*x3).*(x1 - x2).*(x2 - x3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).^2.*(x2 - x3)./ar.^3) + 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x3).*(y1 - y2)./ar + ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar).*(y1 - y2)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x2 - x3).*(y1 - y2)./ar + ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 + x2 - 2.*x3).*(y1 - y2)./ar - ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x2 - x3).*(y1 - y2)./ar + ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar).*(y1 - y3)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 + x2 - 2.*x3).*(y1 - y2)./ar - ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - x3).*(y1 - y2)./ar + ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar).*(y2 - y3)./ar.^2 + 4.*(y1 - y2).*((xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy)./ar - 1)./ar + 4.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y2)./ar.^2);	
    Mb113=-1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x1 - x2).*(y1 + y2 - 2.*y3)./ar + ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x1 - x2).*(y1 + y2 - 2.*y3)./ar + ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x1 - x2).*(y1 - y3)./ar - ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x1 - x2).*(y1 - y3)./ar - ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x1 - x2).*(y2 - y3)./ar - ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x1 - x2).*(y2 - y3)./ar - ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 + 4.*(x1 - x2).*((xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy)./ar - 1)./ar + 4.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x2)./ar.^2) - 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 + y2 - 2.*y3).*(y1 - y2).*(y1 - y3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).^2.*(y1 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 + y2 - 2.*y3).*(y1 - y2).*(y2 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).^2.*(y2 - y3)./ar.^3);	
    Mb114=-1./2.*Db.*(y1 - y2)./ar;	
    Mb115=1./2.*Db.*nu.*(x1 - x2)./ar;	
    Mb21=-1./2.*Db.*nu.*(y2 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar - 1./2.*Db.*(x2 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar;	
    Mb22=1./8.*Db.*nu.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).*((2.*x1 - x2 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).*((2.*x1 - x2 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3).*((x1 - x2).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).*((x1 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 - x2).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y2 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y2 - y3)./ar.^2 + 4.*(y2 - y3).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy)./ar - 1)./ar + 4.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y2 - y3)./ar.^2) + 1./8.*Db.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(2.*x1 - x2 - x3).*(x1 - x2).*(x2 - x3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(2.*x1 - x2 - x3).*(x1 - x3).*(x2 - x3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(x2 - x3).^2./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*(x2 - x3).^2./ar.^3);	
    Mb23=1./8.*Db.*nu.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(2.*y1 - y2 - y3).*(y1 - y2).*(y2 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(2.*y1 - y2 - y3).*(y1 - y3).*(y2 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).*(y2 - y3).^2./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).*(y2 - y3).^2./ar.^3) + 1./8.*Db.*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x2 - x3).*(2.*y1 - y2 - y3)./ar - ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x2 - x3).*(2.*y1 - y2 - y3)./ar - ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x2 - x3).*(y1 - y2)./ar - ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x2 - x3).*(y1 - y2)./ar - ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x2 - x3).*(y1 - y3)./ar - ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x2 - x3).*(y1 - y3)./ar - ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - 4.*(x2 - x3).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy)./ar - 1)./ar - 4.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x2 - x3)./ar.^2);	
    Mb24=-1./2.*Db.*nu.*(y2 - y3)./ar;	
    Mb25=1./2.*Db.*(x2 - x3)./ar;	
    Mb26=1./2.*Db.*nu.*(y1 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar + 1./2.*Db.*(x1 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar;	
    Mb27=-1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - 2.*x2 + x3).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y1 - y2)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x2 - x3).*(y1 - y3)./ar - ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x2).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y1 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x2 - x3).*(y1 - y3)./ar - ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar).*(y1 - y3)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 - x2).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - 2.*x2 + x3).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y2 - y3)./ar.^2 - 4.*(y1 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy)./ar + 1)./ar - 4.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y3)./ar.^2) - 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(x1 - 2.*x2 + x3).*(x1 - x3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(x1 - x3).^2./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - 2.*x2 + x3).*(x1 - x3).*(x2 - x3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).^2.*(x2 - x3)./ar.^3);	
    Mb28=-1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).*(y1 - 2.*y2 + y3).*(y1 - y3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).*(y1 - y3).^2./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - 2.*y2 + y3).*(y1 - y3).*(y2 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).^2.*(y2 - y3)./ar.^3) + 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x1 - x3).*(y1 - y2)./ar + ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x1 - x3).*(y1 - y2)./ar + ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x1 - x3).*(y1 - 2.*y2 + y3)./ar + ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x1 - x3).*(y1 - 2.*y2 + y3)./ar + ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x1 - x3).*(y2 - y3)./ar + ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x1 - x3).*(y2 - y3)./ar + ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 - 4.*(x1 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy)./ar + 1)./ar - 4.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x3)./ar.^2);	
    Mb29=1./2.*Db.*nu.*(y1 - y3)./ar;	
    Mb210=-1./2.*Db.*(x1 - x3)./ar;	
    Mb211=-1./2.*Db.*nu.*(y1 - y2).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar - 1./2.*Db.*(x1 - x2).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar;	
    Mb212=1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x3).*(y1 - y2)./ar + ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar).*(y1 - y2)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x2 - x3).*(y1 - y2)./ar + ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 + x2 - 2.*x3).*(y1 - y2)./ar - ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x2 - x3).*(y1 - y2)./ar + ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar).*(y1 - y3)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 + x2 - 2.*x3).*(y1 - y2)./ar - ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - x3).*(y1 - y2)./ar + ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar).*(y2 - y3)./ar.^2 + 4.*(y1 - y2).*((xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy)./ar - 1)./ar + 4.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y2)./ar.^2) - 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 + x2 - 2.*x3).*(x1 - x2).*(x1 - x3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).^2.*(x1 - x3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 + x2 - 2.*x3).*(x1 - x2).*(x2 - x3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).^2.*(x2 - x3)./ar.^3);	
    Mb213=-1./8.*Db.*nu.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 + y2 - 2.*y3).*(y1 - y2).*(y1 - y3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2).^2.*(y1 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 + y2 - 2.*y3).*(y1 - y2).*(y2 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).^2.*(y2 - y3)./ar.^3) - 1./8.*Db.*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x1 - x2).*(y1 + y2 - 2.*y3)./ar + ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x1 - x2).*(y1 + y2 - 2.*y3)./ar + ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x1 - x2).*(y1 - y3)./ar - ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x1 - x2).*(y1 - y3)./ar - ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x1 - x2).*(y2 - y3)./ar - ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x1 - x2).*(y2 - y3)./ar - ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 + 4.*(x1 - x2).*((xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy)./ar - 1)./ar + 4.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x2)./ar.^2);	
    Mb214=-1./2.*Db.*nu.*(y1 - y2)./ar;	
    Mb215=1./2.*Db.*(x1 - x2)./ar;	
    Mb31=-1./4.*Db.*(nu - 1).*((y2 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar + (x2 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar);	
    Mb32=1./16.*Db.*(nu - 1).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(2.*x1 - x2 - x3).*(x2 - x3).*(y1 - y2)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*(x2 - x3).*(y1 - y2)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(2.*x1 - x2 - x3).*(x2 - x3).*(y1 - y3)./ar.^3 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(x2 - x3).*(y1 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((2.*x1 - x2 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((2.*x1 - x2 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x1 - x2).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x1 - x2).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x1 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x1 - x3).*(y2 - y3)./ar + ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(x2 - x3).*(y2 - y3)./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*(x2 - x3).*(y2 - y3)./ar.^3 + 4.*(x2 - x3).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy)./ar - 1)./ar + 4.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x2 - x3)./ar.^2);	
    Mb33=1./16.*Db.*(nu - 1).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x2 - x3).*(2.*y1 - y2 - y3)./ar - ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x2 - x3).*(y1 - y3)./ar - ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x2 - x3).*(2.*y1 - y2 - y3)./ar - ((x2 - x3).*(y2 + y3) - (x2 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x2 - x3).*(y1 - y2)./ar - ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y1 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x2 - x3).*(y1 - y2)./ar - ((x2 - x3).*(y1 + y2) - (x1 + x2).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y2 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x2 - x3).*(y1 - y3)./ar - ((x2 - x3).*(y1 + y3) - (x1 + x3).*(y2 - y3) + 2.*x3.*y2 - 2.*x2.*y3)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(2.*y1 - y2 - y3).*(y2 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*(2.*y1 - y2 - y3).*(y2 - y3)./ar.^3 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*(y1 - y2).*(y2 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*(y1 - y2).*(y2 - y3)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(y1 - y3).*(y2 - y3)./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*(y1 - y3).*(y2 - y3)./ar.^3 - 4.*(y2 - y3).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy)./ar - 1)./ar - 4.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y2 - y3)./ar.^2);	
    Mb34=-1./4.*Db.*(nu - 1).*(x2 - x3)./ar;	
    Mb35=1./4.*Db.*(nu - 1).*(y2 - y3)./ar;	
    Mb36=1./4.*Db.*(nu - 1).*((y1 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar + (x1 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar);	
    Mb37=1./16.*Db.*(nu - 1).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x1 - x2).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x1 - x2).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x1 - 2.*x2 + x3).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x1 - 2.*x2 + x3).*(y1 - y3)./ar - ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x2 - x3).*(y1 - y3)./ar - ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x2 - x3).*(y1 - y3)./ar - ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - 2.*x2 + x3).*(x1 - x3).*(y1 - y2)./ar.^3 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x3).*(x2 - x3).*(y1 - y2)./ar.^3 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(x1 - x3).*(y1 - y3)./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*(x2 - x3).*(y1 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(x1 - x3).*(y2 - y3)./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - 2.*x2 + x3).*(x1 - x3).*(y2 - y3)./ar.^3 + 4.*(x1 - x3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy)./ar + 1)./ar + 4.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x3)./ar.^2);	
    Mb38=-1./16.*Db.*(nu - 1).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x3).*(y1 - 2.*y2 + y3)./ar + ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x3).*(y1 - y2)./ar + ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y1 - y3)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*(y1 - y2).*(y1 - y3)./ar.^3 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*(y1 - y2).*(y1 - y3)./ar.^3 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(y1 - 2.*y2 + y3).*(y1 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*(y1 - 2.*y2 + y3).*(y1 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).*((x1 - x3).*(y2 - y3)./ar + ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).*((x1 - x3).*(y2 - y3)./ar + ((x2 + x3).*(y1 - y3) - 2.*x3.*y1 - (x1 - x3).*(y2 + y3) + 2.*x1.*y3)./ar)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 - x3).*(y1 - y2)./ar + ((x1 - x3).*(y1 + y2) - (x1 + x2).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - x3).*(y1 - 2.*y2 + y3)./ar + ((x1 - x3).*(y1 + y3) - (x1 + x3).*(y1 - y3) + 2.*x3.*y1 - 2.*x1.*y3)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(y1 - y3).*(y2 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*(y1 - y3).*(y2 - y3)./ar.^3 + 4.*(y1 - y3).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy)./ar + 1)./ar + 4.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y3)./ar.^2);	
    Mb39=1./4.*Db.*(nu - 1).*(x1 - x3)./ar;	
    Mb310=-1./4.*Db.*(nu - 1).*(y1 - y3)./ar;	
    Mb311=-1./4.*Db.*(nu - 1).*((y1 - y2).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3)./ar.^2)./ar + (x1 - x2).*((xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y2)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(y1 - y3)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y2 - y3)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y2 - y3)./ar.^2)./ar);	
    Mb312=-1./16.*Db.*(nu - 1).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*((x1 + x2 - 2.*x3).*(y1 - y2)./ar - ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*((x1 + x2 - 2.*x3).*(y1 - y2)./ar - ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*((x1 - x3).*(y1 - y2)./ar + ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*((x1 - x3).*(y1 - y2)./ar + ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*((x2 - x3).*(y1 - y2)./ar + ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*((x2 - x3).*(y1 - y2)./ar + ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(x1 - x3).*(y1 - y2)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(x2 - x3).*(y1 - y2)./ar.^3 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 + x2 - 2.*x3).*(x1 - x2).*(y1 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x2).*(x2 - x3).*(y1 - y3)./ar.^3 + (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 + x2 - 2.*x3).*(x1 - x2).*(y2 - y3)./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x2).*(x1 - x3).*(y2 - y3)./ar.^3 - 4.*(x1 - x2).*((xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy)./ar - 1)./ar - 4.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x2)./ar.^2);	
    Mb313=1./16.*Db.*(nu - 1).*((xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x2).*(y1 - y3)./ar - ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar).*(y1 - y2)./ar.^2 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x3).*(y1 + y2 - 2.*y3).*(y1 - y2)./ar.^3 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x2 - x3).*(y1 + y2 - 2.*y3).*(y1 - y2)./ar.^3 - (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*((x1 - x2).*(y1 + y2 - 2.*y3)./ar + ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar).*(y1 - y3)./ar.^2 + (xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*(x1 - x2).*(y1 - y2).*(y1 - y3)./ar.^3 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x2 - x3).*(y1 - y2).*(y1 - y3)./ar.^3 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(y1 - y2).*((x1 - x2).*(y2 - y3)./ar - ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y3).*((x1 - x2).*(y2 - y3)./ar - ((x2 + x3).*(y1 - y2) - 2.*x2.*y1 - (x1 - x2).*(y2 + y3) + 2.*x1.*y2)./ar)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*((x1 - x2).*(y1 + y2 - 2.*y3)./ar + ((x1 - x2).*(y1 + y2) - (x1 + x2).*(y1 - y2) + 2.*x2.*y1 - 2.*x1.*y2)./ar).*(y2 - y3)./ar.^2 + (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*((x1 - x2).*(y1 - y3)./ar - ((x1 + x3).*(y1 - y2) - (x1 - x2).*(y1 + y3) - 2.*x2.*y1 + 2.*x1.*y2)./ar).*(y2 - y3)./ar.^2 - (xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*(x1 - x2).*(y1 - y2).*(y2 - y3)./ar.^3 - (xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(x1 - x3).*(y1 - y2).*(y2 - y3)./ar.^3 - 4.*(y1 - y2).*((xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy)./ar - 1)./ar - 4.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*(y1 - y2)./ar.^2);	
    Mb314=-1./4.*Db.*(nu - 1).*(x1 - x2)./ar;	
    Mb315=1./4.*Db.*(nu - 1).*(y1 - y2)./ar;
   
    %Elements of Shear matrix
    Ts11=0;	
    Ts12=0;	
    Ts13=0;	
    Ts14=1./2.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*Ds./ar;	
    Ts15=0;	
    Ts16=0;	
    Ts17=0;	
    Ts18=0;	
    Ts19=-1./2.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*Ds./ar;	
    Ts110=0;	
    Ts111=0;	
    Ts112=0;	
    Ts113=0;	
    Ts114=1./2.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*Ds./ar;	
    Ts115=0;	
    Ts21=0;	
    Ts22=0;	
    Ts23=0;	
    Ts24=0;	
    Ts25=1./2.*(xx.*(y2 - y3) - x3.*y2 + x2.*y3 - (x2 - x3).*yy).*Ds./ar;	
    Ts26=0;	
    Ts27=0;	
    Ts28=0;	
    Ts29=0;	
    Ts210=-1./2.*(xx.*(y1 - y3) - x3.*y1 + x1.*y3 - (x1 - x3).*yy).*Ds./ar;	
    Ts211=0;	
    Ts212=0;	
    Ts213=0;	
    Ts214=0;	
    Ts215=1./2.*(xx.*(y1 - y2) - x2.*y1 + x1.*y2 - (x1 - x2).*yy).*Ds./ar;    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Original implementation from ZHUANG paper
%     % vetor deslocamento
%     U1 = U(it1)'; %transversal displacement node 1 
%     U2 = U(it1+np)'; %rotation in x node 1
%     U3 = U(it1+2*np)'; %rotation in y node 1
%     U4 = U(it1+3*np)'; %transversal shear strain in x node 1
%     U5 = U(it1+4*np)'; %transversal shear strain in y node 1
%     U6 = U(it2)'; %transversal displacement node 2
%     U7 = U(it2+np)'; %rotation in x node 2
%     U8 = U(it2+2*np)'; %rotation in y node 2
%     U9 = U(it2+3*np)'; %transversal shear strain in x node 2
%     U10 = U(it2+4*np)'; %transversal shear strain in y node 2
%     U11 = U(it3)'; %transversal displacement node 3
%     U12 = U(it3+np)'; %rotation in x node 3
%     U13 = U(it3+2*np)'; %rotation in y node 3
%     U14 = U(it3+3*np)'; %transversal shear strain in x node 3
%     U15 = U(it3+4*np)'; %transversal shear strain in y node 3

    % vetor deslocamento
    U1 = U(it1)'; %transversal displacement node 1 
    U3 = U(it1+np)'; %rotation in x node 1
    U2 = U(it1+2*np)'; %rotation in y node 1
    U5 = U(it1+3*np)'; %transversal shear strain in x node 1
    U4 = U(it1+4*np)'; %transversal shear strain in y node 1
    U6 = U(it2)'; %transversal displacement node 2
    U8 = U(it2+np)'; %rotation in x node 2
    U7 = U(it2+2*np)'; %rotation in y node 2
    U10 = U(it2+3*np)'; %transversal shear strain in x node 2
    U9 = U(it2+4*np)'; %transversal shear strain in y node 2
    U11 = U(it3)'; %transversal displacement node 3
    U13 = U(it3+np)'; %rotation in x node 3
    U12 = U(it3+2*np)'; %rotation in y node 3
    U15 = U(it3+3*np)'; %transversal shear strain in x node 3
    U14 = U(it3+4*np)'; %transversal shear strain in y node 3
    
    
    % Global s = [Mxx, Myy, Mxy, Sx, Sy]
    s=zeros(5,nt);
    s(1,:) = Mb11.*U1+Mb12.*U2+Mb13.*U3+Mb14.*U4+Mb15.*U5+Mb16.*U6+Mb17.*U7+Mb18.*U8+Mb19.*U9+Mb110.*U10+Mb111.*U11+Mb112.*U12+Mb113.*U13+Mb114.*U14+Mb115.*U15;
    s(2,:) = Mb21.*U1+Mb22.*U2+Mb23.*U3+Mb24.*U4+Mb25.*U5+Mb26.*U6+Mb27.*U7+Mb28.*U8+Mb29.*U9+Mb210.*U10+Mb211.*U11+Mb212.*U12+Mb213.*U13+Mb214.*U14+Mb215.*U15;
    s(3,:) = Mb31.*U1+Mb32.*U2+Mb33.*U3+Mb34.*U4+Mb35.*U5+Mb36.*U6+Mb37.*U7+Mb38.*U8+Mb39.*U9+Mb310.*U10+Mb311.*U11+Mb312.*U12+Mb313.*U13+Mb314.*U14+Mb315.*U15;
    s(4,:) = Ts11.*U1+Ts12.*U2+Ts13.*U3+Ts14.*U4+Ts15.*U5+Ts16.*U6+Ts17.*U7+Ts18.*U8+Ts19.*U9+Ts110.*U10+Ts111.*U11+Ts112.*U12+Ts113.*U13+Ts114.*U14+Ts115.*U15;
    s(5,:) = Ts21.*U1+Ts22.*U2+Ts23.*U3+Ts24.*U4+Ts25.*U5+Ts26.*U6+Ts27.*U7+Ts28.*U8+Ts29.*U9+Ts210.*U10+Ts211.*U11+Ts212.*U12+Ts213.*U13+Ts214.*U14+Ts215.*U15;

end
end