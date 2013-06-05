%
% FDTD 1D in free space with CPML 
%
% AUTHOR:skiloop
% EMAIL:skiloop@126.com
%

clc
close all;

%% some constants
mu_0 = 1.2566370614359173e-06;   %
eps_0= 8.8541878176203892e-12;   %
C=299792458.0;      % speed of light

%% wave definition
amptidute=100;
frequency=110E9;
T=1/frequency;
lambda=C*T;
t0=4*T;
omega=2*pi*frequency;

%% FDTD variables
if exist('zZoneSize','var')==0
    zZoneSize=4;
else if zZoneSize<0
    zZoneSize=4;
    end
end
if exist('tZoneSize','var')==0
    tZoneSize=4;
else if zZoneSize<0
    tZoneSize=4;
    end
end
domainLength=zZoneSize*lambda;
totalTime=tZoneSize*T;
numberCellsPerWavelength=200;

%% CPML parameters
pmlWidth=40;
pmlOrder=4;
epsR=1;
sigmaMax=1;
kappaMax=15;
alphaMax=0.24;
alphaOrder=1;

%% domain definition
dz=lambda/numberCellsPerWavelength;
dt=dz/2/C;
totalTimeStep=floor(totalTime/dt);
nz=floor(domainLength/dz+2*pmlWidth);
nzp1=nz+1;
nzm1=nz-1;
ksource = floor(nzp1/2);
%% FDTD EM field arrarys
Ex=zeros(1,nzp1);
Hy=zeros(1,nz);

%% Coeficient for EM field updating
Cexe=1;
Cexhy=-dt/eps_0/dz;
Chyh=1;
Chyex=-dt/mu_0/dz;

%% CPML arrays
Psi_Ezy_zn=zeros(1,pmlWidth);
Psi_Ezy_zp=zeros(1,pmlWidth);
Psi_Hzx_zn=zeros(1,pmlWidth);
Psi_Hzx_zp=zeros(1,pmlWidth);

%% initial PML update coeficients
sigmaOpt=sigmaMax*(pmlOrder+1)/(sqrt(epsR)*150*pi*dz);
% zn side
rho_e=((pmlWidth:-1:1)-0.75)/pmlWidth;
rho_m=((pmlWidth:-1:1)-0.25)/pmlWidth;
sigma_e=sigmaOpt*abs(rho_e).^pmlOrder;
sigma_m=sigmaOpt*abs(rho_m).^pmlOrder;
kappa_e=1+(kappaMax-1)*abs(rho_e).^pmlOrder;
kappa_m=1+(kappaMax-1)*abs(rho_m).^pmlOrder;
alpha_e=alphaMax*abs(rho_e).^pmlOrder;
alpha_m=alphaMax*abs(rho_m).^pmlOrder;
cpml_b_e_n=exp((-dt/eps_0)...
    *((sigma_e./kappa_e)+alpha_e));
cpml_a_e_n=1/dz*(cpml_b_e_n-1).*sigma_e ...
    ./(kappa_e.*(sigma_e+kappa_e.*alpha_e));
cpml_b_m_n=exp((-dt/eps_0)...
    *((sigma_m./kappa_m)+alpha_m));
cpml_a_m_n=1/dz*(cpml_b_m_n-1).*sigma_m ...
    ./(kappa_m.*(sigma_m+kappa_m.*alpha_m));
b=0;%sigma_e*dt/2/eps_0;
Cexe_zn=(1-b)./(1+b);
CPsi_Ezy_zn=-dt./(1+b)/eps_0;
Cexhy_zn=CPsi_Ezy_zn./kappa_e/dz;
b=0;%sigma_m*dt/2/mu_0;
Chyh_zn=(1-b)./(1+b);
CPsi_Hzx_zn=-dt./(1+b)/mu_0;
Chyex_zn=CPsi_Hzx_zn./kappa_m/dz;
% zp side
rho_e=((1:1:pmlWidth)-0.75)/pmlWidth;
rho_m=((1:1:pmlWidth)-0.25)/pmlWidth;
sigma_e=sigmaOpt*abs(rho_e).^pmlOrder;
sigma_m=sigmaOpt*abs(rho_m).^pmlOrder;
kappa_e=1+(kappaMax-1)*abs(rho_e).^pmlOrder;
kappa_m=1+(kappaMax-1)*abs(rho_m).^pmlOrder;
alpha_e=alphaMax*abs(rho_e).^pmlOrder;
alpha_m=alphaMax*abs(rho_m).^pmlOrder;
cpml_b_e_p=exp((-dt/eps_0)...
    *((sigma_e./kappa_e)+alpha_e));
cpml_a_e_p=1/dz*(cpml_b_e_p-1).*sigma_e ...
    ./(kappa_e.*(sigma_e+kappa_e.*alpha_e));
cpml_b_m_p=exp((-dt/eps_0)...
    *((sigma_m./kappa_m)+alpha_m));
cpml_a_m_p=1/dz*(cpml_b_m_p-1).*sigma_m ...
    ./(kappa_m.*(sigma_m+kappa_m.*alpha_m));
b=0;%sigma_e*dt/2/eps_0;
Cexe_zp=(1-b)./(1+b);
CPsi_Ezy_zp=-dt./(1+b)/eps_0;
Cexhy_zp=CPsi_Ezy_zp./kappa_e/dz;
b=0;%sigma_m*dt/2/mu_0;
Chyh_zp=(1-b)./(1+b);
CPsi_Hzx_zp=-dt./(1+b)/mu_0;
Chyex_zp=CPsi_Hzx_zp./kappa_m/dz;

%% FDTD loop

% some constants
%====== Gaussian Source ===========
dtDivEps0DivDz=dt/eps_0/dz;
muSource=dtDivEps0DivDz*amptidute * -2.0 /T/T;
%====== Sine Source ===========
% dtDivEps0DivDz=dt/eps_0/dz;
% muSource=dtDivEps0DivDz*amptidute * 2.0*pi*omega;

% initial plot
figure;
h=plot(Ex);
%set(gca,'ylim',[-2 2.5]*1e15);
set(gca,'xlim',[pmlWidth (nz-pmlWidth)]);
grid on;

% point to test performace of CPML
ic=ksource+numberCellsPerWavelength;
cEx=zeros(1,totalTimeStep);

% fdtd loop
for n=1:totalTimeStep
    %============================
    %update Hy
    %============================
    Psi_Hzx_zn=cpml_b_m_n.*Psi_Hzx_zn+cpml_a_m_n.*(Ex(2:pmlWidth+1)-Ex(1:pmlWidth));
    Psi_Hzx_zp=cpml_b_m_p.*Psi_Hzx_zp+cpml_a_m_p.*(Ex(nzp1-pmlWidth+1:nzp1)-Ex(nzp1-pmlWidth:nz));
    if sum(isnan(Psi_Hzx_zp))>0
        display('nan found');
    end
    Hy(1:pmlWidth)=Chyh_zn.*Hy(1:pmlWidth)+Chyex_zn.*(Ex(2:pmlWidth+1)-Ex(1:pmlWidth))+CPsi_Hzx_zn.*Psi_Hzx_zn;
    if sum(isnan(Hy))>0
        display('nan found');
    end
    Hy(nz-pmlWidth+1:nz)=Chyh_zp.*Hy(nz-pmlWidth+1:nz)+Chyex_zp.*(Ex(nzp1-pmlWidth+1:nzp1)-Ex(nzp1-pmlWidth:nz))+CPsi_Hzx_zp.*Psi_Hzx_zp;
    if sum(isnan(Hy))>0
        display('nan found');
    end
    % non pml region
    Hy(pmlWidth+1:nz-pmlWidth)=Chyh* Hy(pmlWidth+1:nz-pmlWidth)+Chyex*( Ex(pmlWidth+2:nzp1-pmlWidth)-Ex(pmlWidth+1:nz-pmlWidth)); 
    
    %===========================
    % update Ex
    %===========================
    Psi_Ezy_zn=cpml_b_e_n.*Psi_Ezy_zn+cpml_a_e_n.*(Hy(2:pmlWidth+1)-Hy(1:pmlWidth));
    Psi_Ezy_zp=cpml_b_e_p.*Psi_Ezy_zp+cpml_a_e_p.*(Hy(nzp1-pmlWidth:nz)-Hy(nz-pmlWidth:nzm1));
    if sum(isnan(Psi_Ezy_zp))>0
        display('nan found');
    end
    Ex(2:pmlWidth+1)=Cexe_zn.*Ex(2:pmlWidth+1)+Cexhy_zn.*(Hy(2:pmlWidth+1)-Hy(1:pmlWidth))+CPsi_Ezy_zn.*Psi_Ezy_zn;
    
    if sum(isnan(Ex))>0
        display('nan found');
    end
    Ex(nzp1-pmlWidth:nz)=Cexe_zp.*Ex(nzp1-pmlWidth:nz)+Cexhy_zp.*(Hy(nzp1-pmlWidth:nz)-Hy(nz-pmlWidth:nzm1))+CPsi_Ezy_zp.*Psi_Ezy_zp;
    if sum(isnan(Ex))>0
        display('nan found');
    end
    % non pml region
    Ex(pmlWidth+2:nz-pmlWidth)=Cexe* Ex(pmlWidth+2:nz-pmlWidth)+Cexhy*( Hy(pmlWidth+2:nz-pmlWidth)-Hy(pmlWidth+1:nzm1-pmlWidth)); 
    
    %==========================
    % update source
    %==========================
    Ex(ksource)=Ex(ksource)+muSource *(n * dt - t0)...
                    * exp(-(((n * dt - t0) / T).^2)); % Differentiated Gaussian pulse
%     Ex(ksource)=Ex(ksource)-amptidute *sin((n * dt - t0) * omega)/dz; % Sine source
                
    %==========================
    % update figure
    %==========================
    if mod(n,30)==0||n==totalTimeStep
        set(h,'YData',Ex);
        title(gca,strcat('time step :' ,int2str(n)));
        pause(0.2);
    end
    
    %=========================
    % sample watched field
    %=========================
    cEx(n)=Ex(ic);    
end