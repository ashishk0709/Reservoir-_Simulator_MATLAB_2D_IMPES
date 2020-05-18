% function PostProcessing
load('P_all.mat');
load('SW_all.mat');
load('Jw_all.mat');
load('Jo_all.mat');
load('QWLw_all.mat');
load('QWLo_all.mat');
load('QWL_all.mat');
load('Pwfw_all.mat');
load('Pwfo_all.mat');
load('QWLt_all.mat');

Pwfw_all(:,1)=[];
Pwfw_all(:,1825:1840)=[];

QWLw_all(:,1:3)=[];
QWLo_all(:,1:3)=[];
QWLt_all(:,1:3)=[];
QWL_all(:,1:3)=[];

Nx=numerical.Nx;
Ny=numerical.Ny;
Tf=numerical.Tf;
wL=wells.L;
nw=wells.nw;
nhzw=size(wells.hzG);
nhzw=nhzw(1);
RW=[0 500 1826 3987];      % time (days) to be plotted
dx=reservoir.dx;
dy=reservoir.dy;
x=dx/2:dx:reservoir.Lx-dx/2;
y=dy/2:dy:reservoir.Ly-dy/2;
t1=linspace(0,Tf,Tf+1);

t1(:,1)=[];
t11=t1(:,3:end);
t3=t11;
t1(:,1825:1840)=[];
t2=t1;
t1=linspace(0,Tf,Tf+1);

%... Pressure Plot at three times
fig=0;
for i=1:length(RW)
    t=RW(i);
    P=P_all(:,t+1);
    P(P==0)=NaN;
    P=(reshape(P,[Nx,Ny]))';
    Pmax=max(max(P));
    Pmin=min(min(P));
    fig=fig+1;
    h1=figure(fig);
    surf(x,y,P)
    zlim([Pmin,Pmax]);
    title1=['Pressure at time=',num2str(t),' Days'];
    title(title1,'FontSize',18);
    xlabel({'X distance (ft)'},'FontSize',14);
    ylabel({'Y distance(ft)'},'FontSize',14);
    view(0,90);
    colorbar;
    for j=1:nw
        x1=wells.X(j);
        y1=wells.Y(j);
        disp1=['WELL-',num2str(j)];
        text(x1,y1,0,disp1);
    end
    title1=['01-',title1];
    print(h1,'-dpng','-r600',title1);
end


%... Saturation Plot at three times
for i=1:length(RW)
    t=RW(i);
    sw1=SW_all(:,t+1);
    sw1(sw1==0)=NaN;
    sw1=(reshape(sw1,[Nx,Ny]))';
    sw1max=max(max(sw1));
    sw1min=min(min(sw1));
    fig=fig+1;
    h1=figure(fig);
    surf(x,y,sw1)
    zlim([sw1min,sw1max]);
    title1=['Saturation at time=',num2str(t),' Days'];
    title(title1,'FontSize',18);
    xlabel({'X distance (ft)'},'FontSize',14);
    ylabel({'Y distance(ft)'},'FontSize',14);
    view(0,90);
    colorbar;
    for j=1:nw
        x1=wells.X(j);
        y1=wells.Y(j);
        disp1=['WELL-',num2str(j)];
        text(x1,y1,0,disp1);
    end
    title1=['02-',title1];
    print(h1,'-dpng','-r600',title1);
end

CM =hsv(nw);  
lin={'-','--',':','-','--',':','-','--'};

%... Plot for well bottom-hole pressures for each well
fig=fig+1;
h1=figure(fig);
le=[];
for i=1:nw    
    bhp1=Pwfw_all(i,:);
    plot(t2,bhp1,'color',CM(i,:),'LineWidth',2.5,'linestyle',lin{i});
    le1={['Well',num2str(i),' BHP']};
    le=[le
        le1];
    hold on;
end
xlabel({'Time (days)'},'FontSize',14);
ylabel({'Pressure (psi)'},'FontSize',14);
title1=['Bottomhole Pressure'];
title(title1,'FontSize',18);
legend(le);
title1=['03-',title1];
pause
print(h1,'-dpng','-r600',title1);


% plot for rate

fig=fig+1;
h1=figure(fig);
le=[];
for i=1:nw
    qw1=abs(QWL_all(i,:));
    plot(t3,qw1,'color',CM(i,:),'LineWidth',2.5,'linestyle',lin{i});
    le1={['Well',num2str(i),' Flow rate']};
    le=[le
        le1];
    hold on;
end
xlabel({'Time (days)'},'FontSize',14);
ylabel({'Total Flow rate (bbl/day)'},'FontSize',14);
ylim([0,2000]);
legend(le);
title1=['Total Liquid Flow rates for each well'];
title(title1,'FontSize',18);
legend(le);
title1=['06-',title1];
print(h1,'-dpng','-r600',title1);

fig=fig+1;
h1=figure(fig);
le=[];
for i=1:nw
    qw1=abs(QWLo_all(i,:));
    plot(t3,qw1,'color',CM(i,:),'LineWidth',2.5,'linestyle',lin{i});
    le1={['Well',num2str(i),' Flow rate']};
    le=[le
        le1];
    hold on;
end
xlabel({'Time (days)'},'FontSize',14);
ylabel({'Oil Flow rate (bbl/day)'},'FontSize',14);
ylim([0,2000]);
legend(le);
title1=['Oil Flow rates for each well'];
title(title1,'FontSize',18);
legend(le);
title1=['07-',title1];
print(h1,'-dpng','-r600',title1);