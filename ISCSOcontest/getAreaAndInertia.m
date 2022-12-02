function rg = getAreaAndInertia()
%     syms r t A I
    close all
    r = sym("r","positive");
    t = sym("t","positive");
    A = sym("A","positive");
    I = sym("I","positive");
    re = r + t;
    ri = r - t;
    eq(1) = A - pi*(re^2 - ri^2);
    eq(2) = I - 1/4*pi*(re^4-ri^4);
    rg = solve(eq,[r,t]);
    rs = rg.r;
    rs = matlabFunction(rs);
    ts = rg.t;
    ts = matlabFunction(ts);
    [Si,indexOpt,h,indexAll] = DataSections();
    Av = Si(:,1);
    Iv = Si(:,2);
    rv = rs(Av,Iv);
    tv = ts(Av,Iv);
    isFeasible = rv - tv > 0;
    rf = rv(isFeasible);
    tf = tv(isFeasible);
    Af = [Av;Av];
    Af = Af(isFeasible);
    If = [Iv;Iv];
    If = If(isFeasible);
    rd = linspace(min(rf),max(rf),20);
    td = linspace(min(tf),max(tf),20);
    [rg,tg] = meshgrid(rd,td);
    re = rg + tg;
    ri = rg - tg;
    A = @(re,ri) pi*(re.^2 - ri.^2);
    I = @(re,ri) 1/4*pi*(re.^4-ri.^4);
    Ag = A(re,ri);
    Ig = I(re,ri);
    ref = rf + tf;
    rif = rf - tf; 
    r2 = (re.^2 + ri.^2)./2;
    r2f = (ref.^2 + rif.^2)./4;
    r2 = linspace(min(r2f(:)),max(r2f(:)),20);
    Ag2 = linspace(min(Af(:)),max(Af(:)),20);
    [r2,Ag2] = meshgrid(r2,Ag2);
    Ig2 = Ag2.*r2;
    Aef = pi*ref.^2;
    Aif = pi*rif.^2;
    [Ae,Ai] = meshgrid(Aef,Aif);
    Aei = Ae - Ai;
    Iei = Aei.*(Ae + Ai)./(4*pi);
    h = h./5;
    plotWithRadAndThick(rg,tg,Ag,Ig,rf,tf,Af,If,indexOpt,h,'r','t')
    plotWithRadAndThick(r2,Ag2,Ag2,Ig2,r2f,Af,Af,If,indexOpt,h,'r^2','Area')
%     plotWithRadAndThick(Ae,Ai,Aei,Iei,Aef,Aif,Af,If,indexOpt,h,'Ae','Ai')
    index = 1:length(Af);
    notIndexOpt = setdiff(index,indexOpt);

    figure()
    hold on
    plot(Af(notIndexOpt),If(notIndexOpt),'ro','MarkerFaceColor', 'r','MarkerSize',3)
    for i = 1:length(h)
        plot(Af(indexOpt(i)),If(indexOpt(i)),'ko','MarkerFaceColor', 'k','MarkerSize',h(i))
    end

    xlabel('Area')
    ylabel('Inertia')
    hold off

    %% Comparing A/I-length values
    % for our solution and the optimal one
%     close all
    load("l.mat");
    load("sect.mat");
    AI = Si(indexAll,1)./Si(indexAll,2);
    figure()
    plot(AI,l,'ro','MarkerFaceColor', 'r')
    xlabel('A/I')
    ylabel('Bar length')

    
    figure()
    Area = Si(sect,1);
    Inertia = Si(sect,2);
    
    subplot(1,2,1)
    hold on
    plot(Area,Inertia,'ko','MarkerFaceColor', 'k')
    plot(Si(:,1),Si(:,2),'ro','MarkerFaceColor', 'r','MarkerSize',3)
    for i = 1:length(h)
        plot(Area(indexOpt(i)),Inertia(indexOpt(i)),'ko','MarkerFaceColor', 'k','MarkerSize',h(i))
        plot(Si(indexAll(i),1),Si(indexAll(i),2),'ro','MarkerFaceColor', 'r','MarkerSize',h(i))
    end
    hold off
    xlabel('Area')
    ylabel('Inertia')
    title('El nostre')

    subplot(1,2,2)
    hold on
    plot(Si(indexAll,1),Si(indexAll,2),'ko','MarkerFaceColor', 'k') %optim
    plot(Si(:,1),Si(:,2),'ro','MarkerFaceColor', 'r','MarkerSize',3) %nostre
    for i = 1:length(h)
        plot(Area(indexOpt(i)),Inertia(indexOpt(i)),'ko','MarkerFaceColor', 'k','MarkerSize',h(i))
        plot(Si(indexAll(i),1),Si(indexAll(i),2),'ro','MarkerFaceColor', 'r','MarkerSize',h(i))
    end
    hold off
    xlabel('Area')
    ylabel('Inertia')
    title('Òptim')

%     close all
    figure()
    hold on
    plot(Area,Inertia,'ro')
    plot(Si(indexAll,1),Si(indexAll,2),'b*')
    hold on
    for i = 1:length(h)
        plot(Area(indexOpt(i)),Inertia(indexOpt(i)),'ro','MarkerSize',h(i))
        plot(Si(indexAll(i),1),Si(indexAll(i),2),'b*','MarkerSize',h(i))
    end
%     plot(Si(:,1),Si(:,2),'ok')
    hold off

end

function plotWithRadAndThick(rg,tg,Ag,Ig,rf,tf,Af,If,indexOpt,h,xlabeln,ylabeln)
    index = 1:length(Af);
    notIndexOpt = setdiff(index,indexOpt);
    figure()
    subplot(1,2,1)
    hold on
    surf(rg,tg,Ag,'FaceAlpha',0.2)
    plot3(rf(notIndexOpt),tf(notIndexOpt),Af(notIndexOpt),'ro','MarkerFaceColor', 'r','MarkerSize',3)
    for i = 1:length(h)
        plot3(rf(indexOpt(i)),tf(indexOpt(i)),Af(indexOpt(i)),'ko','MarkerFaceColor', 'k','MarkerSize',h(i))
    end
    ylabel(ylabeln)
    xlabel(xlabeln)

    title('Area')

    subplot(1,2,2)
    hold on
    surf(rg,tg,Ig,'FaceAlpha',0.2)
    plot3(rf(notIndexOpt),tf(notIndexOpt),If(notIndexOpt),'ro','MarkerFaceColor', 'r','MarkerSize',3)
    for i = 1:length(h)
        plot3(rf(indexOpt(i)),tf(indexOpt(i)),If(indexOpt(i)),'ko','MarkerFaceColor', 'k','MarkerSize',h(i))
    end
    ylabel(ylabeln)
    xlabel(xlabeln)
    title('Inertia')
    hold off

end