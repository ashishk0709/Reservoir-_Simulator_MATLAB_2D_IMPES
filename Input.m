function [reservoir,numerical,fluid,wells]=Input

%% .. Buckley Leverett ..
reservoir.Lx=6000;                               % Reservoir length (ft) in X direction
reservoir.Ly=7500;                               % Reservoir length (ft) in Y direction

reservoir.depth=importdata('PJ1-depth.txt');     % depth of grid block
reservoir.Kx=importdata('PJ1-Permeability.txt'); % Permeability in X direction
reservoir.phi=importdata('PJ1-Porosity.txt');    % Porosity
reservoir.th=importdata('PJ1-Thickness.txt');    % Thickness
[Ny,Nx]=size(reservoir.depth);                   % Number of grid block

reservoir.th=reshape(reservoir.th',[1,Nx*Ny]);
reservoir.Kx=reshape(reservoir.Kx',[1,Nx*Ny]);
reservoir.depth=reshape(reservoir.depth',[1,Nx*Ny]);
reservoir.phi=reshape(reservoir.phi',[1,Nx*Ny]);

reservoir.KR=0.15;                                  % Ky/Kx (Perm ratio)
reservoir.Ky=reservoir.KR*reservoir.Kx;          % Permeability in Y direction
reservoir.Kz=reservoir.Kx;                       % Permeability in Z direction
reservoir.dx=(reservoir.Lx/Nx);                  % Grid size in X direction
reservoir.dy=(reservoir.Ly/Ny);                  % Grid size in Y direction

reservoir.dWOC=-7474.45;
reservoir.PwWOC=4500;


%.... GRAVITY & CAPILLARY PRESSURE calculation selection 
reservoir.gc=-6.939E-3;  % lb/ft3 density to psi/ft conversion factor-- Gravity Switch (Make this to zero)
reservoir.PcS=1;                   % Switch for capillary pressure inclusion (ZERO=> Pc OFF)

%% .......... FLUID PROPERTIES .........
fluid.cr=1E-6;
fluid.cw=3.215E-6;
fluid.co=0E-6;
fluid.muw=0.383211;
fluid.muo=1.03;
fluid.Bw=1.02812;
fluid.Bo=1;
fluid.dw=62.4;
fluid.do=53;

%%
% %Fluid Data vs. Pressure
fluid.Pref=[14.7;47.05;79.4;111.76;144.11;176.46;208.82;241.17;273.53;305.88;338.23;370.59;402.94;435.29;467.65;500;1400;2300;3200;4100;5000;6000;8000;10000;12000;14000]; %Reference Pressure for Bo,Co and Muo
fluid.muoRef=[2.34;2.31;2.28;2.25;2.22;2.19;2.16;2.12;2.09;2.06;2.03;2;1.97;1.95;1.92;1.89;1.33;1.03;0.84;0.72;0.63;0.58;0.51;0.45;0.4;0.35]; 
fluid.BoRef=[1.05;1.05;1.05;1.05;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.09;1.11;1.14;1.18;1.21;1.243;1.31;1.37;1.44;1.51]; 
fluid.coRef=[3.0e-05;3.0e-05;3.0e-05;3.00e-05;3.00e-05;3.00e-05;3.00e-05;3.00e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.00e-05;3.00e-05;3.00e-05;3.0e-05;2.66e-05;1.39e-05;9.08e-06;6.58e-06;5.09e-06;4.6e-06;4.0e-06;3.4e-06;2.8e-06;2.2e-06]; 


%% .......... NUMERICAL PROPERTIES .........
numerical.Nx=Nx;       % (i) Number of grid in X direction
numerical.Ny=Ny;       % (j) Number of grid in Y direction
numerical.dt=1;        % Time step size 
numerical.Tf=3987;      % STOP Time

%% 
%....... WELL PROPERTIES..........

wells.nw=6;                   % number of wells
wells.IP=[-1;-1;-1;-1;0;0];         % 1 for injector -1 for producer 0 for constt BHP

wells.type=[1;1;2;2;1;1];       %Well_Type (Vertical=1,horizontal-X=2,horizontal-Y=3)
wells.constraint=[1;1;1;1;2;2]; %(constant rate=1,constant BHP=2) 
wells.X=[3675;3825;2550;2250;1125;450];      % X co-ordinate
wells.Y=[5600;3600;4400;2700;1100;3100];      % Y co-ordinate
wells.rate=[-300;-150;-400;-250;0;0];  % wells rate constraints
wells.BHP=[0;0;0;0;100;150];     % wells BHP constraints
wells.rw=[0.25;0.25;0.25;0.25;0.25;0.25];     % wellbore radius
wells.hzL=[NaN;NaN;225;225;NaN;NaN];    % If horizontal well then length
wells.S=[0;0;0;0;0;0];  % well skin factor


%% 
% BOUNDARY CONDITIONS
     
reservoir.BCL=0;  % Boundary condition-left ( 0 for Neumann & 1 for Drichlet)
reservoir.BCR=0;  % Boundary condition-right (0 for Neumann & 1 for Drichlet)
reservoir.BCT=0;  % Boundary condition-top ( 0 for Neumann & 1 for Drichlet)
reservoir.BCB=0;  % Boundary condition-bottom (0 for Neumann & 1 for Drichlet)

reservoir.PBL=0;  % Pressure at Boundary LEFT
reservoir.PBR=0;  % Pressure at Boundary RIGHT
reservoir.PBT=0;  % Pressure at Boundary TOP
reservoir.PBB=0;  % Pressure at Boundary BOTTOM

reservoir.QL=0;    % Flow rate at Boundary LEFT
reservoir.QR=0;    % Flow rate at Boundary RIGHT
reservoir.QT=0;    % Flow rate at Boundary TOP
reservoir.QB=0;    % Flow rate at Boundary BOTTOM