function [swi,Po]=sw_init

[reservoir,numerical,fluid,~]=Input;
Nx=numerical.Nx;
Ny=numerical.Ny;

phoW=fluid.dw;                                       % Density of Water in lb/ft3
phoO=fluid.do;                                         % Density of Oil in lb/ft3
gc=-6.939E-3;                                    % lb/ft3 to psi/ft conversion factor
phogW=phoW*gc;
phogO=phoO*gc;
dWOC=reservoir.dWOC;                                   % Depth of water oil contact
PwWOC=reservoir.PwWOC;                                      % Pressure of water at WOC
swWOC=1;                                         % Sw at Water Oil Contact
reservoir.depth(reservoir.depth==0)=NaN;

Pe=3.5;
swcon=0.2;
lambda=2;
PcWOC=Pe*((swWOC-swcon)/(1-swcon))^(-1/lambda);
PoWOC=PwWOC+PcWOC;

Pw=PwWOC+phogW*(reservoir.depth-dWOC);
Po=PoWOC+phogO*(reservoir.depth-dWOC);
Pc=Po-Pw;
swi=((1-swcon)*((Pc./Pe).^(-lambda)))+swcon;
swi(swi>=1)=1;
soi=1-swi;

% fig=0;
% %...FOR PLOTTING
% dx=reservoir.dx;
% dy=reservoir.dy;
% x=dx/2:dx:reservoir.Lx-dx/2;
% y=dy/2:dy:reservoir.Ly-dy/2;
% 
% fig=fig+1;
% h1=figure(fig);
% soi=(reshape(soi,[Nx,Ny]))';
% surf(x,y,soi)
% zlim([0,1]);
% title1=['Initial Oil Saturation'];
% title(title1,'FontSize',18);
% xlabel({'X distance (ft)'},'FontSize',14);
% ylabel({'Y distance(ft)'},'FontSize',14);
% view(0,90);
% colorbar;
% % print(h1,'-dpng','-r600',title1);
% 
% fig=fig+1;
% h1=figure(fig);
% Po=(reshape(Po,[Nx,Ny]))';
% surf(x,y,Po)
% title1=['Oil Phase Pressure'];
% title(title1,'FontSize',18);
% xlabel({'X distance (ft)'},'FontSize',14);
% ylabel({'Y distance(ft)'},'FontSize',14);
% view(0,90);
% colorbar;
% % print(h1,'-dpng','-r600',title1);
% 
% fig=fig+1;
% h1=figure(fig);
% Pw=(reshape(Pw,[Nx,Ny]))';
% surf(x,y,Pw)
% title1=['Water Phase Pressure'];
% title(title1,'FontSize',18);
% xlabel({'X distance (ft)'},'FontSize',14);
% ylabel({'Y distance(ft)'},'FontSize',14);
% view(0,90);
% colorbar;
% % print(h1,'-dpng','-r600',title1);