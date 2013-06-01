%
% FDTD 1D in free space with CPML 
%
% AUTHOR:skiloop
% EMAIL:skiloop@126.com
%

clc
clear;
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

%% FDTD variables
domainLength=10*lambda;
totalTime=20*T;
numberCellsPerWavelength=100;

%% CPML parameters
pmlWidth=10;
pmlOrder=4;
epsR=1;
kappaMax=8;
alphaMax=0.05;
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
Ce1=1;
Ce2=-dt/eps_0/dz;
Ch1=1;
Ch2=-dt/mu_0/dz;

%% CPML arrays
Psi_Ezy_left=zeros(1,pmlWidth);
Psi_Ezy_right=zeros(1,pmlWidth);
Psi_Hzx_left=zeros(1,pmlWidth);
Psi_Hzx_right=zeros(1,pmlWidth);

%% initial PML update coeficients
sigmaMax=(pmlOrder+1)/(sqrt(epsR)*150*pi*dz);
% left side
rho_e=((pmlWidth:-1:1)-0.75)/pmlWidth;
rho_m=((pmlWidth:-1:1)-0.25)/pmlWidth;
sigma_e=sigmaMax*abs(rho_e).^pmlOrder;
sigma_m=(mu_0/eps_0)*sigmaMax*abs(rho_m).^pmlOrder;
kappa_e=1+(kappaMax-1)*abs(rho_e).^pmlOrder;
kappa_m=1+(kappaMax-1)*abs(rho_m).^pmlOrder;
alpha_e=alphaMax*abs(rho_e).^pmlOrder;
alpha_m=(mu_0/eps_0)*alphaMax*abs(rho_m).^pmlOrder;
cpml_b_e_n=exp((-dt/eps_0)...
    *((sigma_e./kappa_e)+alpha_e));
cpml_a_e_n=1/dz*(cpml_b_e_n-1).*sigma_e ...
    .*(kappa_e.*(sigma_e+kappa_e.*alpha_e));
cpml_b_m_n=exp((-dt/mu_0)...
    *((sigma_m./kappa_m)+alpha_m));
cpml_a_m_n=1/dz*(cpml_b_m_n-1).*sigma_m ...
    .*(kappa_m.*(sigma_m+kappa_m.*alpha_m));
b=  sigma_e*dt/2/eps_0;
ca_left=(1-b)./(1+b);
cc_left=-dt./(1+b)/eps_0;
cb_left=cc_left./kappa_e/dz;
b=  sigma_m*dt/2/mu_0;
c1_left=(1-b)./(1+b);
c3_left=-dt./(1+b)/mu_0;
c2_left=c3_left./kappa_m/dz;
% right side
rho_e=((1:1:pmlWidth)-0.75)/pmlWidth;
rho_m=((1:1:pmlWidth)-0.25)/pmlWidth;
sigma_e=sigmaMax*abs(rho_e).^pmlOrder;
sigma_m=(mu_0/eps_0)*sigmaMax*abs(rho_m).^pmlOrder;
kappa_e=1+(kappaMax-1)*abs(rho_e).^pmlOrder;
kappa_m=1+(kappaMax-1)*abs(rho_m).^pmlOrder;
alpha_e=alphaMax*abs(rho_e).^pmlOrder;
alpha_m=(mu_0/eps_0)*alphaMax*abs(rho_m).^pmlOrder;
cpml_b_e_p=exp((-dt/eps_0)...
    *((sigma_e./kappa_e)+alpha_e));
cpml_a_e_p=1/dz*(cpml_b_e_n-1).*sigma_e ...
    .*(kappa_e.*(sigma_e+kappa_e.*alpha_e));
cpml_b_m_p=exp((-dt/mu_0)...
    *((sigma_m./kappa_m)+alpha_m));
cpml_a_m_p=1/dz*(cpml_b_m_p-1).*sigma_m ...
    .*(kappa_m.*(sigma_m+kappa_m.*alpha_m));
b=  sigma_e*dt/2/eps_0;
ca_right=(1-b)./(1+b);
cc_right=-dt./(1+b)/eps_0;
cb_right=cc_right./kappa_e/dz;
b=  sigma_m*dt/2/mu_0;
c1_right=(1-b)./(1+b);
c3_right=-dt./(1+b)/mu_0;
c2_right=c3_right./kappa_m/dz;

%% FDTD loop

% some constants
dtDivEps0DivDz=dt/eps_0/dz;
muSource=dtDivEps0DivDz*amptidute * -2.0 / T/T;

% initial plot
figure;
h=plot(Ex);
set(gca,'ylim',[-2.5 3]*1e15);
grid on;

% point to test performace of CPML
ic=ksource+floor((nzp1-ksource)/2);
cEx=zeros(1,totalTime);

% fdtd loop
for n=1:totalTimeStep
    %============================
    %update Hy
    %============================
    Psi_Hzx_left=cpml_b_m_n.*Psi_Hzx_left+cpml_b_m_n.*(Ex(2:pmlWidth+1)-Ex(1:pmlWidth));
    Psi_Hzx_right=cpml_b_m_p.*Psi_Hzx_right+cpml_b_m_p.*(Ex(nzp1-pmlWidth+1:nzp1)-Ex(nzp1-pmlWidth:nz));
    Hy(1:pmlWidth)=c1_left.*Hy(1:pmlWidth)+c2_left.*(Ex(2:pmlWidth+1)-Ex(1:pmlWidth))+c3_left.*Psi_Hzx_left;
    Hy(nz-pmlWidth+1:nz)=c1_right.*Hy(nz-pmlWidth+1:nz)+c2_right.*(Ex(nzp1-pmlWidth+1:nzp1)-Ex(nzp1-pmlWidth:nz))+c3_right.*Psi_Hzx_right;
    % non pml region
    Hy(pmlWidth+1:nz-pmlWidth)=Ch1* Hy(pmlWidth+1:nz-pmlWidth)+Ch2*( Ex(pmlWidth+2:nzp1-pmlWidth)-Ex(pmlWidth+1:nz-pmlWidth)); 
    
    %===========================
    % update Ex
    %===========================
    Psi_Ezy_left=cpml_b_e_n.*Psi_Ezy_left+cpml_b_e_n.*(Hy(2:pmlWidth+1)-Hy(1:pmlWidth));
    Psi_Ezy_right=cpml_b_e_p.*Psi_Ezy_right+cpml_b_e_p.*(Hy(nzp1-pmlWidth:nz)-Hy(nz-pmlWidth:nzm1));
    Ex(2:pmlWidth+1)=ca_left.*Ex(2:pmlWidth+1)+cb_left.*(Hy(2:pmlWidth+1)-Hy(1:pmlWidth))+cc_left.*Psi_Ezy_left;
    Ex(nzp1-pmlWidth:nz)=ca_right.*Ex(nzp1-pmlWidth:nz)+cb_right.*(Hy(nzp1-pmlWidth:nz)-Hy(nz-pmlWidth:nzm1))+cc_right.*Psi_Ezy_right;
    % non pml region
    Ex(pmlWidth+2:nzp1-pmlWidth)=Ce1* Ex(pmlWidth+2:nzp1-pmlWidth)+Ce2*( Hy(pmlWidth+2:nzp1-pmlWidth)-Hy(pmlWidth+1:nz-pmlWidth)); 
    
    %==========================
    % update source
    %==========================
    Ex(ksource)=Ex(ksource)+muSource* ((n * dt - t0) ) ...
                    * exp(-(((n * dt - t0) / T).^2)); % Differentiated Gaussian pulse
                
    %==========================
    % update figure
    %==========================
    if mod(n,10)==0
        set(h,'YData',Ex);
        title(gca,strcat('time step :' ,int2str(n)));
        pause(0.2);
    end
    
    %=========================
    % sample watched field
    %=========================
    cEx(n)=Ex(ic);
    
end

%% display performance of CPML
% caculate analysis solve
t=(1:totalTimeStep)*dt;
delay=(ic-ksource)*dz/C;
aEx=muSource* ((t - t0-delay) ).*exp(-(((t - t0-delay) / T).^2));
relativeError=abs(cEx-aEx)/max(abs(aEx));
dbError=20*log10(relativeError);

figure('NumberTitle','OFF','Name','Time Macthing');
plot(t/1e-9,aEx,t/1e-9,cEx,'LineWidth',2);
xlabel('time (ns)');
ylabel('Ez');
grid on;
title('time domain compare');
legend('analysis','cpml');

figure('NumberTitle','OFF','Name','Relative Error');
semilogy(t/1e-9,relativeError,'LineWidth',2);
xlabel('time (ns)');
ylabel('Relative Error');
grid on;
title('Relative Error');
    
figure('NumberTitle','OFF','Name','DB Error');
plot(t/1e-9,dbError,'LineWidth',2);
xlabel('time (ns)');
ylabel('Error (DB)');
grid on;
title('DB Error');



