clear all
clc
load(strcat('PaperFiguresAndResults/LargePlastic_HighForcing_1.mat'),'bh');
b=copyBroyden(bh{17});
load(strcat('PaperFiguresAndResults/LargePlastic_MediumForcing_2.mat'),'bh');
b=copyBroyden(bh{30});

% nice but small!
load(strcat('PaperFiguresAndResults/LargePlastic_LowForcing_1.mat'),'bh');
b=copyBroyden(bh{13});

% nice!
load(strcat('PaperFiguresAndResults/NarrowPlastic_MediumForcing_2.mat'),'bh');
b=copyBroyden(bh{22});

load(strcat('PaperFiguresAndResults/OnTopPlastic_MediumForcing_2.mat'),'bh');
b=copyBroyden(bh{32});

load(strcat('PaperFiguresAndResults/OnTopPlastic_MediumForcing_4.mat'),'bh');
b=copyBroyden(bh{28});
b=copyBroyden(bh{2});

%% Copy quantities from a certain freq



st=b.record_sensors;
stH=b.record_sensors_harm;


lag=st(9,1);
Res=b.record_residual(:,end);
% Res = b.Zglobal*[stH(:,5);stH(:,6)] + b.Rtarget + b.S2N*[stH(:,1);stH(:,2)] 
% XN = [stH(:,5);stH(:,6)]
% XP = - b.Zglobal\(b.S2N*[stH(:,1);stH(:,2)] + b.Rtarget) = [stH(:,5);stH(:,6)]-b.Zglobal\(Res)
% FP = [stH(:,1);stH(:,2)] 
% FN = b.S2N\(b.Zglobal*[stH(:,5);stH(:,6)] + b.Rtarget) = b.S2N\(Res-FP)

% Zglobal = (obj.T2N'*obj.Znt*obj.T2N)*obj.L2T
% (obj.T2N'*obj.Znt*obj.T2N)*obj.L2T*XP = 
%     b.S2N*[stH(:,1);stH(:,2)]+b.Rtarget
% S2N only lag compensation

% [forS; forC; momS; momC]
FP = (b.T2N')\(b.S2N*[stH(:,1);stH(:,2)]);
FN = (b.T2N')\(b.Zglobal*[stH(:,5);stH(:,6)] + b.Rtarget) ;
% [disS; disC; rotS; rotC]
XN=b.T2N*b.L2T*[stH(:,5);stH(:,6)];
XP=-b.T2N*b.L2T*(b.Zglobal\(b.S2N*[stH(:,1);stH(:,2)]+b.Rtarget));

% Plot Results

figure()
subplot(2,2,1)
% b.OfflinePlotsTime(st(1,:),lag,'Measurement');hold on;
b.OfflinePlotsFourier(-FN(1:end/2),'Virtual Interface');
b.OfflinePlotsFourier(FP(1:end/2),'Physical Interface');
ylabel('[N]');xlabel('nondimensional time');title('Force');%legend();hold off

subplot(2,2,2)
% b.OfflinePlotsTime(st(2,:),lag,'Measurement');hold on;
b.OfflinePlotsFourier(-FN(end/2+1:end),'Virtual Interface');
b.OfflinePlotsFourier(FP(end/2+1:end),'Physical Interface');
ylabel('[Nm]');xlabel('nondimensional time');title('Moment');legend();hold off

subplot(2,2,3)
b.OfflinePlotsFourier(XN(1:end/2),'Virtual Interface');
b.OfflinePlotsFourier(XP(1:end/2),'Physical Interface');
% b.OfflinePlotsTime(st(5,:),lag,'Measurement');hold on;
ylabel('[mm]');xlabel('nondimensional time');title('Displacement');%legend();hold off

subplot(2,2,4)
b.OfflinePlotsFourier((XN(end/2+1:end)),'Virtual Interface');
b.OfflinePlotsFourier((XP(end/2+1:end)),'Physical Interface');
% b.OfflinePlotsTime(st(6,:),lag,'Measurement');hold on;
ylabel('[^o]');xlabel('nondimensional time');title('Rotation');legend();hold off

% subplot(2,3,[3,6])
% semilogy(vecnorm(b.record_residual,2,1),'Marker','s');%xlabel('iterations');ylabel('err');


%%
N_periods = 1;
sampling_rate = round(1/b.time_step);
Fourier_periods = [1:floor(sampling_rate/b.freq*N_periods)]/floor(sampling_rate/b.freq);
iDFT = [sin(2*pi*(1:b.Hmax).*Fourier_periods'),...
        cos(2*pi*(1:b.Hmax).*Fourier_periods')];
figure()
subplot(2,2,2)
plot(iDFT*XN(end/2+1:end),iDFT*XP(end/2+1:end))
xlabel('Rotation at NS [^o]');
ylabel('Rotation at PS [^o]');
subplot(2,2,1)
plot(iDFT*XN(1:end/2),iDFT*XP(1:end/2))
xlabel('Displacement at NS [mm]');
ylabel('Displacement at PS [mm]');
subplot(2,2,4)
plot(-iDFT*FN(end/2+1:end),iDFT*FP(end/2+1:end))
xlabel('Moment at NS [Nmm]');
ylabel('Moment at PS [Nmm]');
subplot(2,2,3)
plot(-iDFT*FN(1:end/2),iDFT*FP(1:end/2))
xlabel('Force at NS [N]');
ylabel('Force at PS [N]');
%%
nH=1;
figure(4)
for iT=1:1
    subplot(2,nH,iT)
    b.OfflinePlotsFourier(b.record_residual(1:2*b.Hmax,1),' ');
    hold on
    b.OfflinePlotsFourier(b.record_residual(1+2*b.Hmax:4*b.Hmax,1),' ');
    xlabel('period');ylabel('initial residual');
    title(strcat('Harmonics up to'," ",num2str(iT)))
    
    subplot(2,nH,nH+iT)
    b.OfflinePlotsFourier(b.record_residual(1:2*b.Hmax,end),' ');
    hold on
    b.OfflinePlotsFourier(b.record_residual(1+2*b.Hmax:4*b.Hmax,end),' ');
    xlabel('period');ylabel('final residual');
end

%
figure(5)
for iT=1:nH
    subplot(2,nH,iT)
    b.OfflinePlotsFourier(b.record_xsol(1:2*b.Hmax,1),' ');
    hold on
    b.OfflinePlotsFourier(b.record_xsol(1+2*b.Hmax:4*b.Hmax,1),' ');
    xlabel('period');ylabel('initial voltage');
    title(strcat('Harmonics up to'," ",num2str(iT)))
    
    subplot(2,nH,nH+iT)
    b.OfflinePlotsFourier(b.record_xsol(1:2*b.Hmax,end),' ');
    hold on
    b.OfflinePlotsFourier(b.record_xsol(1+2*b.Hmax:4*b.Hmax,end),' ');
    xlabel('period');ylabel('final voltage');
end

