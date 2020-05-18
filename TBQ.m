function [T,Tw,D,J,Q,G,d11,d12i,WQW,WQO,wells]=TBQ(t,P,SW,reservoir,fluid,numerical,wells)
format long;

Nx=numerical.Nx;
Ny=numerical.Ny;
depth=reservoir.depth';
gw=fluid.dw*reservoir.gc;           % water gradient in psi/ft
go=fluid.do*reservoir.gc;           % water gradient in psi/ft

Tw=sparse(Nx*Ny,Nx*Ny);
To=sparse(Nx*Ny,Nx*Ny);
D=sparse(Nx*Ny,Nx*Ny);
QR=sparse(Nx*Ny,1);
d11=sparse(Nx*Ny,1);
d12i=sparse(Nx*Ny,1);
d21=sparse(Nx*Ny,1);
d22=sparse(Nx*Ny,1);

[Jw,Jo,WQW,WQO,wells]=wellM(t,P,SW,reservoir,fluid,numerical,wells);
Pc=real(PC(SW,reservoir.PcS));
Pc(isnan(Pc))=0;
Pc(Pc>=1000)=1000;
Pc1=real(PC((SW+0.002),reservoir.PcS));
Pc1(isnan(Pc1))=0;
Pc1(Pc1>=1000)=1000;
dPCdSW=(Pc1-Pc)/0.002;    % Calculation of derivative of Pc w.r.t. Sw
dPCdSW(dPCdSW>1e7)=1e7;
dPCdSW(dPCdSW<-1e7)=-1e7;

for l=1:Nx*Ny
    Vi=reservoir.dx*reservoir.dy*reservoir.th(l);
    LR=mod(l,Nx);   % To check for grid on left/right (if zero then right if ONE then left)
    
    TB=0;           % Grid Not on Bottom
    if(l<=Nx)
    TB=1;           % Grid on Bottom
    end
    TT=0;           % Grid Not on Top
    if(Ny==1)
    TT=1;           % Grid on Top
    elseif(l>=Nx*Ny-Nx+1 && l<=Nx*Ny)
    TT=1;           % Grid on Top
    end
    d11(l,l)=(Vi*SW(l)*reservoir.phi(l)*(fluid.cw+fluid.cr))/(fluid.Bw*numerical.dt);
    d12i(l,l)=1/((Vi*reservoir.phi(l)/(fluid.Bw*numerical.dt))*(1-(SW(l)*reservoir.phi(l)*fluid.cw*dPCdSW(l))));
    d21(l,l)=(Vi*reservoir.phi(l)*(fluid.co(l)+fluid.cr)*(1-SW(l)))/(fluid.Bo(l)*numerical.dt);
    d22(l,l)=(-Vi*reservoir.phi(l))/(fluid.Bo(l)*numerical.dt);
    
    
    if(LR~=1)                                         % Block not on Left
        Tw(l,l-1)=-TXWhalf(l-1,l,P,SW,reservoir,fluid);
        Tw(l,l)=Tw(l,l)-Tw(l,l-1);
        To(l,l-1)=-TXOhalf(l-1,l,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)-To(l,l-1);        
    elseif(LR==1 && reservoir.BCL==1)  % Block on Left Dirichlet BC
        Tw(l,l)=Tw(l,l)+2*TXWhalf(l,l,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)+2*TXOhalf(l,l,P,SW,reservoir,fluid);
        QR(l)=QR(l)+2*TXWhalf(l,l,P,SW,reservoir,fluid)*reservoir.PBL;
    elseif(LR==1 && reservoir.BCL==0)  % Block on Left Neumann BC
        QR(l)=QR(l)+reservoir.QL;
    end
    
    if(LR~=0)                                         % Block not on Right
        Tw(l,l+1)=-TXWhalf(l,l+1,P,SW,reservoir,fluid);
        Tw(l,l)=Tw(l,l)-Tw(l,l+1);
        To(l,l+1)=-TXOhalf(l,l+1,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)-To(l,l+1);
    elseif(LR==0 && reservoir.BCR==1) % Block on Right Dirichlet BC
        Tw(l,l)=Tw(l,l)+2*TXWhalf(l,l,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)+2*TXOhalf(l,l,P,SW,reservoir,fluid);
        QR(l)=QR(l)+2*TXWhalf(l,l,P,SW,reservoir,fluid)*reservoir.PBR;
    elseif(LR==0 && reservoir.BCR==0) % Block on Right Neumann BC
        QR(l)=QR(l)+reservoir.QR;
    end
    
    if(TB~=1)                                         % Block not on Bottom
        Tw(l,l-Nx)=-TYWhalf(l-Nx,l,P,SW,reservoir,fluid);
        Tw(l,l)=Tw(l,l)-Tw(l,l-Nx);
        To(l,l-Nx)=-TYOhalf(l-Nx,l,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)-To(l,l-Nx);        
    elseif(TB==1 && reservoir.BCB==1) % Block on Bottom Dirichlet BC
        Tw(l,l)=Tw(l,l)+2*TYWhalf(l,l,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)+2*TYOhalf(l,l,P,SW,reservoir,fluid);
        QR(l)=QR(l)+2*TYWhalf(l,l,P,SW,reservoir,fluid)*reservoir.PBB;
    elseif(TB==1 && reservoir.BCB==0) % Block on Bottom Neumann BC
        QR(l)=QR(l)+reservoir.QB;
    end
    
    if(TT~=1)                                             % Block not on Top 
        Tw(l,l+Nx)=-TYWhalf(l,l+Nx,P,SW,reservoir,fluid);
        Tw(l,l)=Tw(l,l)-Tw(l,l+Nx);
        To(l,l+Nx)=-TYOhalf(l,l+Nx,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)-To(l,l+Nx);        
    elseif(TT==1 && reservoir.BCT==1)       % Block on Top Dirichlet BC
        Tw(l,l)=Tw(l,l)+2*TYWhalf(l,l,P,SW,reservoir,fluid);
        To(l,l)=To(l,l)+2*TYOhalf(l,l,P,SW,reservoir,fluid);
        QR(l)=QR(l)+2*TYWhalf(l,l,P,SW,reservoir,fluid)*reservoir.PBT;
    elseif(TT==1 && reservoir.BCT==0)       % Block on Top Neumann BC
        QR(l)=QR(l)+reservoir.QT;
    end
    

D(l,l)=((-d22(l,l)*d12i(l,l))*d11(l,l))+d21(l,l);    
end
T=((-d22*d12i)*Tw)+To;
J=((-d22*d12i)*Jw)+Jo;
Q=sparse(((-d22*d12i)*WQW)+WQO -d22*d12i*Jw*Pc);

G=sparse((-d22*d12i*Tw*Pc)+(-d22*d12i*gw*Tw+go*To)*depth);


