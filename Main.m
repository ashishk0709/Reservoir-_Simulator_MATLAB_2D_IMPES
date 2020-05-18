clc
clear
warning off;
format long;

[reservoir,numerical,fluid,wells]=Input;
dt=numerical.dt;         % TIME STEP SIZE in days
Tf=numerical.Tf;         % STOP Time
t=0+dt;
gw=fluid.dw*reservoir.gc;           % water gradient in psi/ft
depth=reservoir.depth';

[SW_old,P_old]=sw_init;
SW_old=SW_old';    % Saturation Initialize 
P_old=P_old';      % Pressure Initialize 

% % Initial Pressure & Saturation for Buckley Leverett 
% P_old=1000*ones(numerical.Nx*numerical.Ny,1);
% SW_old=0.2*ones(numerical.Nx*numerical.Ny,1);

P_old(isnan(P_old))=0;
SW_old(isnan(SW_old))=0;


P_all=[P_old];
SW_all=[SW_old];
Jw_all=zeros(wells.nw,1);
Jo_all=zeros(wells.nw,1);
QWLw_all=zeros(wells.nw,1);
QWLo_all=zeros(wells.nw,1);
Pwfw_all=zeros(wells.nw,1);
Pwfo_all=zeros(wells.nw,1);
QWLt_all=zeros(wells.nw,1);

%....... Calculation for each time step .........
while t<=Tf

    %.... fluid properties updated as per pressure
    fluid.muo=interp1(fluid.Pref,fluid.muoRef,P_old,'cubic','extrap');   %Oil Viscosity cp
    fluid.Bo=interp1(fluid.Pref,fluid.BoRef,P_old,'linear','extrap');    %Oil FVF
    fluid.co=interp1(fluid.Pref,fluid.coRef,P_old,'cubic','extrap');     %Oil Viscosity cp
    
    % T,B,Q Matrix
    [T,Tw,D,J,Q,G,d11,d12i,WQW,WQO,wells]=TBQ(t,P_old,SW_old,reservoir,fluid,numerical,wells);

    % Pressure Solver
    P_new=(T+J+D)\(D*P_old+Q+G);
    P_new(isnan(P_new))=0;
    Pc=real(PC(SW_old,reservoir.PcS));
    Pc(isnan(Pc))=0;
    Pc(Pc>=1000)=1000;
    
    % Saturation equation
    SW_new=SW_old + d12i*( -Tw*(P_new) + Tw*Pc +gw*Tw*depth -d11*(P_new-P_old)+ WQW  );
    SW_new(isnan(SW_new))=0;
    P_old=P_new;
    SW_old=SW_new;
    
    
    % Saving Well Index Values for flow calculation in post-processing
    Jw1=wells.WIW';                     
    Jo1=wells.WIO';
    Jw_all=[Jw_all Jw1];
    Jo_all=[Jo_all Jo1];
    
    %%... Saving Pressure and saturation for each time step
    P_all=[P_all P_new];           % Time and Pressure stored in a matrix for each iteration
    SW_all=[SW_all SW_new];        % Time and Saturation stored in a matrix for each iteration
    display([num2str(floor(100*t/Tf)),' % completed']);   
    t=t+dt;                         % Next Time step
    
    
    %%... Saving BHP and flow rate for each time step    
    
    for i=1:wells.nw
        QWLw(i)=0;              % Water flow rate from particular well
        QWLo(i)=0;              % Oil flow rate from particular well
        J1w=wells.WIW(i);       % Well Index for water
        J1o=wells.WIO(i);       % well Index for Oil
        JT=full(J(wells.L(i),wells.L(i)));
        Pb=P_new(wells.L(i));   % well grid block pressure
        nhzw=numel(find(wells.type==2));           % number of horizontal well
        if (wells.type(i)==1 && wells.constraint(i)==2)         % Vertical well constt BHP
            Pwfw(i)=wells.BHP(i);
            Pwfo(i)=wells.BHP(i);
            QWLw(i)=J1w*(wells.BHP(i)-Pb)/5.615;        % flow rate in bbl/day
            QWLo(i)=J1o*(wells.BHP(i)-Pb)/5.615;
            QWLt(i)=JT*(wells.BHP(i)-Pb)/5.615;
        elseif (wells.type(i)==1 && wells.constraint(i)==1)     % Vertical well constt RATE
            if (wells.IP(i)==1)        % INJECTOR 
                QWLw(i)=wells.rate(i);
                QWLo(i)=0;
                Pwfw(i)=Pb+(5.615*QWLw(i)/J1w);
                Pwfo(i)=0;
                QWLt(i)=wells.rate(i);
            else                        
                qT=wells.rate(i);       % bbl/day
                swW=SW_new(wells.L(i));
                [krw,kro]=KR(swW);
                fw=(krw/fluid.muw)/((krw/fluid.muw)+(kro/fluid.muo(wells.L(i))));            
                QWLw(i)=fw*qT;
                QWLo(i)=(1-fw)*qT;
                Pwfw(i)=Pb+(5.615*QWLw(i)/J1w);
                Pwfo(i)=Pb+(5.615*QWLo(i)/J1o);
                QWLt(i)=wells.rate(i);
            end
        elseif (wells.type(i)==2 && wells.constraint(i)==1)     % Horizontal well constt rate
            if (wells.IP(i)==1)        % INJECTOR 
                QWLw(i)=wells.rate(i);      % bbl/day
                QWLo(i)=0;
                Pwfw(i)=Pb+(5.615*QWLw(i)/J1w/nhzw);
                Pwfo(i)=0;
                QWLt(i)=wells.rate(i);
            else          
                qT=wells.rate(i);       %bbl/day
                swW=SW_new(wells.L(i));
                [krw,kro]=KR(swW);
                fw=(krw/fluid.muw)/((krw/fluid.muw)+(kro/fluid.muo(wells.L(i))));            
                QWLw(i)=fw*qT;      % bbl/day
                QWLo(i)=(1-fw)*qT;  % bbl/day
                Pwfw(i)=Pb+(5.615*QWLw(i)/J1w/nhzw);
                Pwfo(i)=Pb+(5.615*QWLo(i)/J1o/nhzw);
                QWLt(i)=wells.rate(i);
            end
        elseif (wells.type(i)==2 && wells.constraint(i)==2)
            Pwfw(i)=wells.BHP(i);
            Pwfo(i)=wells.BHP(i);
            for j=wells.hzG(1,i):wells.hzG(end,i)
                if(j==0)
                else
                Pb=P_new(j);
                qw1=J1w*(wells.BHP(i)-Pb)/5.615;
                qo1=J1o*(wells.BHP(i)-Pb)/5.615;
                QWLw(i)=QWLw(i)+qw1;
                QWLo(i)=QWLo(i)+qo1;
                end
            end
            QWLt(i)=JT*nhzw*(wells.BHP(i)-Pb)/5.615;
        end
    end
    QWLw_all=[QWLw_all QWLw'];          % Water flow from each well
    QWLo_all=[QWLo_all QWLo'];          % Oil flow from each well
    Pwfw_all=[Pwfw_all Pwfw'];
    Pwfo_all=[Pwfo_all Pwfo'];
    QWLt_all=[QWLt_all QWLt'];
end

P_all(P_all==0)=NaN;
dx=reservoir.dx;
x=dx/2:dx:reservoir.Lx-dx/2;
xD=x/reservoir.Lx;

% %%.. Plot for comparison with Buckley leverett
% t1=[57 114 188];
% tD=(5.615*wells.rate(1)*t1)/(reservoir.th(1)*reservoir.dy*reservoir.Lx*reservoir.phi(3));
% h1=figure(1);
% for i=1:numel(t1)
% t2=t1(i);
% tD1=tD(i);
% [xD2,sw2]=BL(tD1);
% sw1=SW_all(:,t2);
% plot(xD,sw1)
% hold on
% plot(xD2,sw2,'r')
% xlim([0,1]);
% ylim([0,1]);
% xlabel({'xD'},'FontSize',14);
% ylabel({'Sw'},'FontSize',14);
% hold on;
% end
% print(h1,'-dpng','-r600','Buckley');
% hold off;

save('P_all.mat','P_all');      % Save results to a MAT file
save('SW_all.mat','SW_all');      % Save results to a MAT file
save('Jw_all.mat','Jw_all');      % Save results to a MAT file
save('Jo_all.mat','Jo_all');      % Save results to a MAT file
save('Pwfo_all.mat','Pwfo_all');      % Save results to a MAT file
save('Pwfw_all.mat','Pwfw_all');      % Save results to a MAT file
save('QWLw_all.mat','QWLw_all');      % Save results to a MAT file
save('QWLo_all.mat','QWLo_all');      % Save results to a MAT file
QWL_all=QWLw_all+QWLo_all;
save('QWL_all.mat','QWL_all');      % Save results to a MAT file
save('QWLt_all.mat','QWLt_all');      % Save results to a MAT file

PostProcessing