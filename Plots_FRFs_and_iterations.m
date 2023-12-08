clear all
clc
close all
figure()
figure()
figure()
figure()

for iFile=1:2
    load(strcat('PaperFiguresAndResults/LargePlastic_MediumForcing_',num2str(iFile),'.mat'),'bh');
    AFm=[];ABm=[];pFm=[];pBm=[];Omm=[];Fexc=[];N_it=[];
    for iFreq=1:length(bh)
        if bh{iFreq}.Conv
            N_it = [N_it;size(bh{iFreq}.record_residual,2)];
            Lfs=bh{iFreq}.record_sensors_harm(1,5);
            Lfc=bh{iFreq}.record_sensors_harm(2,5);
            Lbs=bh{iFreq}.record_sensors_harm(1,6);
            Lbc=bh{iFreq}.record_sensors_harm(2,6);
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;(Lfs.^2+Lfc.^2).^0.5];
            ABm=[ABm;(Lbs.^2+Lbc.^2).^0.5];
            pFm=[pFm;atan2(Lfs,-Lfc)];
            pBm=[pBm;atan2(Lbs,-Lbc)];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        else
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;0];
            ABm=[ABm;0];
            pFm=[pFm;-pi/2];
            pBm=[pBm;-pi/2];
            N_it = [N_it;100];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        end
    end
    %
    m=1;
    figure(1)
    hold on
    subplot(1,2,1)%subplot(2,2,1)
    hold on
    % plot(Om./(2*pi),AF.*Fexc,'-')
    plot(Omm./(2*pi),m*AFm,'g.')
    % hold off
    title('Front laser amplitude')
    subplot(1,2,2)%subplot(2,2,2)
    hold on
    % plot(Om./(2*pi),pF,'-')
    plot(Omm./(2*pi),pFm,'g.')
    % hold off
    title('Front laser phase')
%     subplot(2,2,3)
%     hold on
%     % plot(Om./(2*pi),AB.*Fexc,'-')
%     plot(Omm./(2*pi),m*ABm,'.')
%     title('Back laser amplitude')
%     subplot(2,2,4)
%     hold on
%     % plot(Om./(2*pi),pB,'-')
%     plot(Omm./(2*pi),pBm,'.')
%     % hold off
%     title('Back laser phase')
    figure(2)
    subplot(1,3,1)
    hold on
    plot(Omm./(2*pi),N_it,'g.')
    ylim([0 100]);xlim([16,20])
    title('Number of iterations')
end


for iFile=1:2
    load(strcat('PaperFiguresAndResults/NarrowPlastic_MediumForcing_',num2str(iFile),'.mat'),'bh');
    AFm=[];ABm=[];pFm=[];pBm=[];Omm=[];Fexc=[];N_it=[];
    for iFreq=1:length(bh)
        if bh{iFreq}.Conv
            N_it = [N_it;size(bh{iFreq}.record_residual,2)];
            Lfs=bh{iFreq}.record_sensors_harm(1,5);
            Lfc=bh{iFreq}.record_sensors_harm(2,5);
            Lbs=bh{iFreq}.record_sensors_harm(1,6);
            Lbc=bh{iFreq}.record_sensors_harm(2,6);
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;(Lfs.^2+Lfc.^2).^0.5];
            ABm=[ABm;(Lbs.^2+Lbc.^2).^0.5];
            pFm=[pFm;atan2(Lfs,-Lfc)];
            pBm=[pBm;atan2(Lbs,-Lbc)];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        else
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;0];
            ABm=[ABm;0];
            pFm=[pFm;-pi/2];
            pBm=[pBm;-pi/2];
            N_it = [N_it;100];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        end
    end
    %
    m=1;
    figure(1)
    hold on
    subplot(1,2,1)%subplot(2,2,1)
    hold on
    % plot(Om./(2*pi),AF.*Fexc,'-')
    plot(Omm./(2*pi),m*AFm,'b.')
    % hold off
    title('Front laser amplitude')
    subplot(1,2,2)%subplot(2,2,2)
    hold on
    % plot(Om./(2*pi),pF,'-')
    plot(Omm./(2*pi),pFm,'b.')
    % hold off
    title('Front laser phase')
%     subplot(2,2,3)
%     hold on
%     % plot(Om./(2*pi),AB.*Fexc,'-')
%     plot(Omm./(2*pi),m*ABm,'.')
%     title('Back laser amplitude')
%     subplot(2,2,4)
%     hold on
%     % plot(Om./(2*pi),pB,'-')
%     plot(Omm./(2*pi),pBm,'.')
%     % hold off
%     title('Back laser phase')
    figure(2)
    hold on
    subplot(1,3,2)
    hold on
    plot(Omm./(2*pi),N_it,'b.')
    ylim([0 100]);xlim([16,20])
    title('Number of iterations')
end


for iFile=1:5
    load(strcat('PaperFiguresAndResults/OnTopPlastic_MediumForcing_',num2str(iFile),'.mat'),'bh');
    AFm=[];ABm=[];pFm=[];pBm=[];Omm=[];Fexc=[];N_it=[];
    for iFreq=1:length(bh)
        if bh{iFreq}.Conv
            N_it = [N_it;size(bh{iFreq}.record_residual,2)];
            Lfs=bh{iFreq}.record_sensors_harm(1,5);
            Lfc=bh{iFreq}.record_sensors_harm(2,5);
            Lbs=bh{iFreq}.record_sensors_harm(1,6);
            Lbc=bh{iFreq}.record_sensors_harm(2,6);
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;(Lfs.^2+Lfc.^2).^0.5];
            ABm=[ABm;(Lbs.^2+Lbc.^2).^0.5];
            pFm=[pFm;atan2(Lfs,-Lfc)];
            pBm=[pBm;atan2(Lbs,-Lbc)];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        else
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;0];
            ABm=[ABm;0];
            pFm=[pFm;-pi/2];
            pBm=[pBm;-pi/2];
            N_it = [N_it;100];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        end
    end
    %
    m=2;
    figure(1)
    hold on
    subplot(1,2,1)%subplot(2,2,1)
    hold on
    % plot(Om./(2*pi),AF.*Fexc,'-')
    plot(Omm./(2*pi),m*AFm,'r.')
    % hold off
    title('Front laser amplitude')
    subplot(1,2,2)%subplot(2,2,2)
    hold on
    % plot(Om./(2*pi),pF,'-')
    plot(Omm./(2*pi),pFm,'r.')
    % hold off
    title('Front laser phase')
%     subplot(2,2,3)
%     hold on
%     % plot(Om./(2*pi),AB.*Fexc,'-')
%     plot(Omm./(2*pi),m*ABm,'.')
%     title('Back laser amplitude')
%     subplot(2,2,4)
%     hold on
%     % plot(Om./(2*pi),pB,'-')
%     plot(Omm./(2*pi),pBm,'.')
%     % hold off
%     title('Back laser phase')
    figure(2)
    hold on
    subplot(1,3,3)
    hold on
    plot(Omm./(2*pi),N_it,'r.')
    ylim([0 100]);xlim([16,20])
    title('Number of iterations')
end


for iFile=1:2
    load(strcat('PaperFiguresAndResults/LargePlastic_LowForcing_',num2str(iFile),'.mat'),'bh');
    AFm=[];ABm=[];pFm=[];pBm=[];Omm=[];Fexc=[];N_it=[];
    for iFreq=1:length(bh)
        if bh{iFreq}.Conv
            N_it = [N_it;size(bh{iFreq}.record_residual,2)];
            Lfs=bh{iFreq}.record_sensors_harm(1,5);
            Lfc=bh{iFreq}.record_sensors_harm(2,5);
            Lbs=bh{iFreq}.record_sensors_harm(1,6);
            Lbc=bh{iFreq}.record_sensors_harm(2,6);
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;(Lfs.^2+Lfc.^2).^0.5];
            ABm=[ABm;(Lbs.^2+Lbc.^2).^0.5];
            pFm=[pFm;atan2(Lfs,-Lfc)];
            pBm=[pBm;atan2(Lbs,-Lbc)];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        else
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;0];
            ABm=[ABm;0];
            pFm=[pFm;-pi/2];
            pBm=[pBm;-pi/2];
            N_it = [N_it;100];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        end
    end
    %
    m=1;
    figure(3)
    hold on
    subplot(1,2,1)%subplot(2,2,1)
    hold on
    % plot(Om./(2*pi),AF.*Fexc,'-')
    plot(Omm./(2*pi),m*AFm,'y.')
    % hold off
    title('Front laser amplitude')
    subplot(1,2,2)%subplot(2,2,2)
    hold on
    % plot(Om./(2*pi),pF,'-')
    plot(Omm./(2*pi),pFm,'y.')
    % hold off
    title('Front laser phase')
%     subplot(2,2,3)
%     hold on
%     % plot(Om./(2*pi),AB.*Fexc,'-')
%     plot(Omm./(2*pi),m*ABm,'.')
%     title('Back laser amplitude')
%     subplot(2,2,4)
%     hold on
%     % plot(Om./(2*pi),pB,'-')
%     plot(Omm./(2*pi),pBm,'.')
%     % hold off
%     title('Back laser phase')
    figure(4)
    subplot(1,3,1)
    hold on
    plot(Omm./(2*pi),N_it,'y.')
    ylim([0 100]);xlim([16,19])
    title('Number of iterations')
end


for iFile=1:2
    load(strcat('PaperFiguresAndResults/LargePlastic_MediumForcing_',num2str(iFile),'.mat'),'bh');
    AFm=[];ABm=[];pFm=[];pBm=[];Omm=[];Fexc=[];N_it=[];
    for iFreq=1:length(bh)
        if bh{iFreq}.Conv
            N_it = [N_it;size(bh{iFreq}.record_residual,2)];
            Lfs=bh{iFreq}.record_sensors_harm(1,5);
            Lfc=bh{iFreq}.record_sensors_harm(2,5);
            Lbs=bh{iFreq}.record_sensors_harm(1,6);
            Lbc=bh{iFreq}.record_sensors_harm(2,6);
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;(Lfs.^2+Lfc.^2).^0.5];
            ABm=[ABm;(Lbs.^2+Lbc.^2).^0.5];
            pFm=[pFm;atan2(Lfs,-Lfc)];
            pBm=[pBm;atan2(Lbs,-Lbc)];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        else
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;0];
            ABm=[ABm;0];
            pFm=[pFm;-pi/2];
            pBm=[pBm;-pi/2];
            N_it = [N_it;100];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        end
    end
    %
    m=1;
    figure(3)
    hold on
    subplot(1,2,1)%subplot(2,2,1)
    hold on
    % plot(Om./(2*pi),AF.*Fexc,'-')
    plot(Omm./(2*pi),m*AFm,'g.')
    % hold off
    title('Front laser amplitude')
    subplot(1,2,2)%subplot(2,2,2)
    hold on
    % plot(Om./(2*pi),pF,'-')
    plot(Omm./(2*pi),pFm,'g.')
    % hold off
    title('Front laser phase')
%     subplot(2,2,3)
%     hold on
%     % plot(Om./(2*pi),AB.*Fexc,'-')
%     plot(Omm./(2*pi),m*ABm,'.')
%     title('Back laser amplitude')
%     subplot(2,2,4)
%     hold on
%     % plot(Om./(2*pi),pB,'-')
%     plot(Omm./(2*pi),pBm,'.')
%     % hold off
%     title('Back laser phase')
    figure(4)
    hold on
    subplot(1,3,2)
    hold on
    plot(Omm./(2*pi),N_it,'g.')
    ylim([0 100]);xlim([16,19])
    title('Number of iterations')
end


for iFile=1:2
    load(strcat('PaperFiguresAndResults/LargePlastic_HighForcing_',num2str(iFile),'.mat'),'bh');
    AFm=[];ABm=[];pFm=[];pBm=[];Omm=[];Fexc=[];N_it=[];
    for iFreq=1:length(bh)
        if bh{iFreq}.Conv
            N_it = [N_it;size(bh{iFreq}.record_residual,2)];
            Lfs=bh{iFreq}.record_sensors_harm(1,5);
            Lfc=bh{iFreq}.record_sensors_harm(2,5);
            Lbs=bh{iFreq}.record_sensors_harm(1,6);
            Lbc=bh{iFreq}.record_sensors_harm(2,6);
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;(Lfs.^2+Lfc.^2).^0.5];
            ABm=[ABm;(Lbs.^2+Lbc.^2).^0.5];
            pFm=[pFm;atan2(Lfs,-Lfc)];
            pBm=[pBm;atan2(Lbs,-Lbc)];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        else
            Fexc=[Fexc;norm(bh{iFreq}.Rtarget)];
            AFm=[AFm;0];
            ABm=[ABm;0];
            pFm=[pFm;-pi/2];
            pBm=[pBm;-pi/2];
            N_it = [N_it;100];
            % disp(strcat(num2str(iFreq)," ",num2str(bh{1}.Conv)))
            Omm=[Omm;bh{iFreq}.freq*2*pi];
        end
    end
    %
    m=1;
    figure(3)
    hold on
    subplot(1,2,1)%subplot(2,2,1)
    hold on
    % plot(Om./(2*pi),AF.*Fexc,'-')
    plot(Omm./(2*pi),m*AFm,'.','Color',[.5,0,1])
    % hold off
    title('Front laser amplitude')
    subplot(1,2,2)%subplot(2,2,2)
    hold on
    % plot(Om./(2*pi),pF,'-')
    plot(Omm./(2*pi),pFm,'.','Color',[.5,0,1])
    % hold off
    title('Front laser phase')
%     subplot(2,2,3)
%     hold on
%     % plot(Om./(2*pi),AB.*Fexc,'-')
%     plot(Omm./(2*pi),m*ABm,'.')
%     title('Back laser amplitude')
%     subplot(2,2,4)
%     hold on
%     % plot(Om./(2*pi),pB,'-')
%     plot(Omm./(2*pi),pBm,'.')
%     % hold off
%     title('Back laser phase')
    figure(4)
    hold on
    subplot(1,3,3)
    hold on
    plot(Omm./(2*pi),N_it,'.','Color',[.5,0,1])
    ylim([0 100]);xlim([16,19])
    title('Number of iterations')
end


