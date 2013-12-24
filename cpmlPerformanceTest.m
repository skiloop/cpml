% 
% Test Performance of cpml
%
clear;

%% reference computing
zZoneSize=10;
tZoneSize=50;
fdtd1dCPML;
aEx=cEx;
save refRes;

% %% no reflect computing
% zZoneSize=100;
% tZoneSize=50;
% fdtd1dCPML;
% save noReflect;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% performance analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=(1:totalTimeStep)*dt;

%%
load refRes aEx;
load noReflect cEx;
relativeError=abs(cEx-aEx)/max(abs(aEx));
dbError=20*log10(relativeError);
% comparision of Ex in time
figure('NumberTitle','OFF','Name','Time Macthing');
plot(t/1e-9,aEx,t/1e-9,cEx,'--','LineWidth',2);
xlabel('time (ns)');
ylabel('Ez');
grid on;
title('time domain compare');
legend('analysis','cpml');
% print -depsc -tiff -r300 cpml1dTimeCompare_r
% print -dtiff -r300 cpml1dTimeCompare_r
%% Relative Error
figure('NumberTitle','OFF','Name','Relative Error');
semilogy(t/1e-9,relativeError,'LineWidth',2);
xlabel('time (ns)');
ylabel('Relative Error');
grid on;
title('Relative Error');
% print -depsc -tiff -r300 cpml1dRelativeError_r
% print -dtiff -r300 cpml1dRelativeError_r
% Error in DB
figure('NumberTitle','OFF','Name','DB Error');
plot(t/1e-9,dbError,'LineWidth',2);
xlabel('time (ns)');
ylabel('Error (DB)');
grid on;
title('DB Error');
% print -depsc -tiff -r300 cpml1dDBError_r
% print -dtiff -r300 cpml1dDBError_r
% show max error
maxDBError=max(dbError)
maxRelError=max(relativeError)
