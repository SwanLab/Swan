function [w,x,y]=triquad(order)
    pospoints=[1;3;4;6;7;12;13;16;19;25];
    points=pospoints(order);
    switch order
    case 1
      wrt=.5;
    case 2
      wrt=[.1;.3;0];
    case 3
      wrt=[.1;.1;.3;0];
    case 4
      wrt=[.05;.1;.2;.4;pi;0];
    case 5
      wrt=[.1;.05;.1;.2;.4;pi;0];
    case 6
      wrt=[.05;.1;.1;.1;.45;.3;.3;.15;0;pi/6;-pi/6;0];
    case 7
      wrt=[-.07;.08;.1;.1;.03;.45;.3;.3;.1;0;pi/6;-pi/6;0];
    case 8
      wrt=[.07;.015;.05;.01;.01;.05;.2;.25;.4;.4;.5;pi;0;pi/14;-pi/14;0];
    case 9
      wrt=[.04;.015;.04;.05;.02;.02;.01;.27;.18;.25;.36;.36;.5;pi;pi;0;pi/8;-pi/8;0];
    case 10
      wrt=[.04;.015;.04;.04;.01;.01;.005;.005;.02;.26;.2;.2;.35;.35;.51;.51;.387;pi;pi/16;-pi/16;pi/8;-pi/8;pi/64;-pi/64;0];
    end

    disp('Solver output:');
    [wrt,R,exr]=fsolve(@(wrt)triwrteq(wrt,order),wrt,optimset('TolFun',1e-16,'TolX',1e-16));

    wxy=wrt2wxy(wrt,order);

    w=2*wxy(1:points);
    x=wxy(points+(1:points));
    y=wxy(2*points+(1:points));

    T=coordtrM([0;1;.5],[0;0;(1-.25)^.5]);

    xy1=T*[wxy(points+(1:points))';wxy(2*points+(1:points))';ones(1,points)];
    edg=T*[0 1 0 0;0 0 1 0;ones(1,4)];
    close all;
    fig=figure('visible','off');
    hold on;axis equal;
    plot(xy1(1,:),xy1(2,:),'o',edg(1,:),edg(2,:));

    cr=unique(floor(sum((xy1(1:2,:)-[.5;3^.5/6]).^2,1).^.5*1000)/1000);
    for k=1:length(cr)
      t=linspace(0,2*pi);
      plot(cr(k)*cos(t)+.5,cr(k)*sin(t)+3^.5/6);
    end
    set(fig,'visible','on');
end
    

    function [R]=triwrteq(wrt,order)
    wxy=wrt2wxy(wrt,order);
    R=triwxeq(wxy,order);
    if order==10
      R=[R;wrt(end-6)+wrt(end-5);wrt(end-4)+wrt(end-3);wrt(end-1)+wrt(end-2);wrt(end);];
    elseif order>=6&&order<=9
      R=[R;wrt(end-1)+wrt(end-2);wrt(end);wrt(end-3)];
    elseif order>1
      R=[R;wrt(end)];
    end
    end
    
    
    function wxy=wrt2wxy(wrt,order)
    pospoints=[1;3;4;6;7;12;13;16;19;25];
    points=pospoints(order);
    wxy=zeros(3*points,1);
    if rem(points,3)==0
      for k=1:points/3
        wxy(3*(k-1)+(1:3))=wrt(k);
      end
      wxy(points+(1:2*points))=rt2points(wrt(points/3+1:end));
    else
      wxy(1)=wrt(1);
      for k=1:(points-1)/3
        wxy(3*(k-1)+1+(1:3))=wrt(k+1);
      end
      wxy([points+1 2*points+1])=1/3;
      if points>1
        wxy([points+1+(1:points-1) 2*(points)+1+(1:points-1)])=rt2points(wrt(1+(points-1)/3+1:end));
      end
    end
    end
    function xy=rt2points(rt)
      lrt=length(rt);
      xy=zeros(3*lrt,1);
      r=rt(1:lrt/2)(:);
      t=rt(lrt/2+(1:lrt/2));
      angs=t+(0:2*pi/3:4*pi/3)+pi/2;
      x=r.*cos(angs)+.5;
      y=r.*sin(angs)+3^.5/6;
      T=coordtrM([0;1;.5],[0;0;(1-.25)^.5]);
      for k=1:lrt/2
        xyt=T\[x(k,:);y(k,:);ones(1,3)];
        xy([3*(k-1)+(1:3) 3*(k-1)+3/2*lrt+(1:3)])=xyt(1:2,1:3)'(:);
      end
    end
    function T=coordtrM(x,y)
    ud=[0;0];
    vd=[1;0];
    wd=[0;1];
    Md=[vd-ud wd-ud ud;0 0 1];
    x=x(:);
    y=y(:);
    M=[x(2:end).'-x(1) x(1);
    y(2:end).'-y(1) y(1);
    zeros(1,2),1];
    T=(M/Md);
    end
    function [R,J]=triwxeq(wxy,order)
    pospoints=[1;3;4;6;7;12;13;16;19;25];
    points=pospoints(order);
    neq=cumsum(1:order+1)(end);
    R=zeros(neq,1);
    J=zeros(neq,3*points);
    intv=zeros(neq,1);
    i=1;
    for o=0:order
      for t=0:o
        xp=o-t;
        yp=t;
        R(i,1)=wxy(1:points).'*(wxy(points+(1:points)).^xp.*wxy(2*points+(1:points)).^yp);
        J(i,:)=[(wxy(points+(1:points)).^xp.*wxy(2*points+(1:points)).^yp).'  (wxy(1:points).*wxy(points+(1:points)).^(max(xp-1,1))*xp.*wxy(2*points+(1:points)).^yp).' (wxy(1:points).*wxy(points+(1:points)).^xp.*wxy(2*points+(1:points)).^(max(yp-1,0))*yp).'];
        intv(i++,1)=sum(bincoeff(yp+1,0:yp+1).*(-1).^((0:yp+1))./((0:yp+1)+xp+1))/(yp+1);
      end
    end
    R-=intv;
end