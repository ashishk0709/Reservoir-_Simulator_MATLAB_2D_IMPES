function [Jw,Jo,WQW,WQO,wells]=wellM(t,P,SW,reservoir,fluid,numerical,wells)
format long;

Nx=numerical.Nx;
Ny=numerical.Ny;
Jw=sparse(Nx*Ny,Nx*Ny);
Jo=sparse(Nx*Ny,Nx*Ny);
WQW=sparse(Nx*Ny,1);
WQO=sparse(Nx*Ny,1);

    
if(t>=1826)
    wells.IP=[0;0;0;0;1;1];         % 1 for injector -1 for producer 0 for constt BHP
    wells.constraint=[2;2;2;2;1;1]; %(constant rate=1,constant BHP=2)
    wells.BHP=[100;75;50;100;0;0];     % wells BHP constraints
    wells.rate=[0;0;0;0;1500;1500];  % wells rate constraints
    wells.S=[0;0;0;0;0;0];  % well skin factor
end

for i=1:wells.nw
    wells.I(i)=ceil(wells.X(i)/reservoir.dx);
    wells.J(i)=ceil(wells.Y(i)/reservoir.dy);
    wells.L(i)=(wells.J(i)-1)*numerical.Nx+wells.I(i);
end

for i=1:wells.nw
    re=0.14*(reservoir.dx^2+reservoir.dy^2)^0.5;
    wg=wells.L(i);
    [krw,kro]=KR(SW(wg));
    kx=reservoir.Kx(wg);
    ky=reservoir.Ky(wg);
    kz=reservoir.Kz(wg);
    dh=reservoir.th(wg);
    rw=wells.rw(i);
    s=wells.S(i);
    if(wells.type(i)==1)
        wells.WIW(i)=(6.33E-3*2*pi*dh*krw*sqrt(kx*ky))/(fluid.muw*fluid.Bw*(log(re/rw)+s));
        wells.WIO(i)=(6.33E-3*2*pi*dh*kro*sqrt(kx*ky))/(fluid.muo(wg)*fluid.Bo(wg)*(log(re/rw)+s));
    elseif(wells.type(i)==2)
        wells.WIW(i)=(6.33E-3*2*pi*krw*reservoir.dx*sqrt(ky*kz))/(fluid.muw*fluid.Bw*(log(re/rw)+s));
        wells.WIO(i)=(6.33E-3*2*pi*kro*reservoir.dx*sqrt(ky*kz))/(fluid.muo(wg)*fluid.Bo(wg)*(log(re/rw)+s));
    end
end

for i=1:wells.nw
    if(wells.type(i)==2)        % Horizontal well in X direction
        hzLf=wells.L(i)+ceil(wells.hzL(i)/reservoir.dx)-1;
        k=1;
        wells.hzG(k,i)=wells.L(i);
        wells.hzTw(k,i)=TXWhalf(wells.L(i),wells.L(i),P,SW,reservoir,fluid);
        wells.hzTo(k,i)=TXOhalf(wells.L(i),wells.L(i),P,SW,reservoir,fluid);
        for j=wells.L(i):hzLf-1
            k=k+1;
            wells.hzG(k,i)=wells.hzG(k-1,i)+1;
            wells.hzTw(k,i)=TXWhalf(wells.hzG(k,i),wells.hzG(k,i),P,SW,reservoir,fluid);
            wells.hzTo(k,i)=TXOhalf(wells.hzG(k,i),wells.hzG(k,i),P,SW,reservoir,fluid);
        end
    end
end

for i=1:wells.nw
    if(wells.type(i)==2)
        p=size(wells.hzG);
        p=p(1);
        for j=1:p
            wells.hzqfw(j,i)=wells.hzTw(j,i)/sum(wells.hzTw(:,i));  % fraction of water flow from each perf
            if(isnan(wells.hzqfw(j,i)))
                wells.hzqfw(j,i)=0;
            end
            wells.hzqfo(j,i)=wells.hzTo(j,i)/sum(wells.hzTo(:,i));  % fraction of oil flow from each perf
            if(isnan(wells.hzqfo(j,i)))
                wells.hzqfo(j,i)=0;
            end
        end
    end
end


%....... WELL MATRIX..........
for i=1:wells.nw
    w_g=wells.L(i);
    [krw,kro]=KR(SW(w_g));
    rw=wells.rw(i);
    s=wells.S(i);  % well skin factor
    kx=reservoir.Kx(w_g);
    ky=reservoir.Ky(w_g);
    dx=reservoir.dx;
    dh=reservoir.th(w_g);
    fw=(krw/fluid.muw)/((krw/fluid.muw)+(kro/fluid.muo(w_g)));
    ip=wells.IP(i);
    if(ip==1)                   % Constant Rate Injector                       
        Wqw=5.615*wells.rate(i);
        Wqo=0;
    else                        % Constant Rate Producer
        Wqw=5.615*fw*wells.rate(i);
        Wqo=5.615*(1-fw)*wells.rate(i);
    end
    
    %.... VERTICAL WELL ....
    if(wells.type(i)==1 && wells.constraint(i)==2)      % Constant BHP
        Pwf=wells.BHP(i);
        re=0.14*(reservoir.dx^2+reservoir.dy^2)^0.5;
        Jw(w_g,w_g)=Jw(w_g,w_g)+(6.33E-3*2*pi*dh*krw*sqrt(kx*ky))/(fluid.muw*fluid.Bw*(log(re/rw)+s));
        Jo(w_g,w_g)=Jo(w_g,w_g)+(6.33E-3*2*pi*dh*kro*sqrt(kx*ky))/(fluid.muo(w_g)*fluid.Bo(w_g)*(log(re/rw)+s));
        WQW(w_g)=WQW(w_g)+Jw(w_g,w_g)*Pwf;
        WQO(w_g)=WQO(w_g)+Jo(w_g,w_g)*Pwf;
    elseif(wells.type(i)==1 && wells.constraint(i)==1)  % Constant Rate 
        WQW(w_g)=WQW(w_g)+Wqw;
        WQO(w_g)=WQO(w_g)+Wqo;
    end
    
    %.... HORIZONTAL WELL in X direction.... 
    if(wells.type(i)==2)
        if(wells.constraint(i)==2)  % Constant BHP
            Pwf=wells.BHP(i);
            for p=wells.hzG(1,i):wells.hzG(end,i)
                if(p==0)
                else
                ky=reservoir.Ky(p);
                kz=reservoir.Kz(p);
                [krw,kro]=KR(SW(p));
                re=0.14*(reservoir.dx^2+reservoir.dy^2)^0.5;
                Jw(p,p)=Jw(p,p)+(6.33E-3*2*pi*reservoir.dx*krw*sqrt(ky*kz))/(fluid.muw*fluid.Bw*(log(re/rw)+s));
                WQW(p)=WQW(p)+Jw(p,p)*Pwf;
                Jo(p,p)=Jo(p,p)+(6.33E-3*2*pi*reservoir.dx*kro*sqrt(ky*kz))/(fluid.muo(p)*fluid.Bo(p)*(log(re/rw)+s));
                WQO(p)=WQO(p)+Jo(p,p)*Pwf;                
                end
            end
        elseif(wells.constraint(i)==1) % Constant Rate
            ip=wells.IP(i);
                if(ip==1)                   % Constant Rate Injector 
                    Wqw=5.615*wells.rate(i);
                    Wqo=0;
                else                        % Constant Rate Producer
                    Wqw=5.615*fw*wells.rate(i);
                    Wqo=5.615*(1-fw)*wells.rate(i);
                end         
            j=1;
            for p=wells.hzG(1,i):wells.hzG(end,i)
                if(p==0)
                else              
                [krw,kro]=KR(SW(p));
                fw=(krw/fluid.muw)/((krw/fluid.muw)+(kro/fluid.muo(p)));
                qgw=Wqw*wells.hzqfw(j,i);  % water flow from the perf                
                WQW(p)=WQW(p)+qgw; 
                qgo=Wqo*wells.hzqfo(j,i);  % oil flow from the perf
                WQO(p)=WQO(p)+qgo;
                j=j+1;
                end
            end
        end            
    end  
end