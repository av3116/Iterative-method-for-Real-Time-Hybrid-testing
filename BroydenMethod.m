classdef BroydenMethod
    %BROYDENMETHOD Iterative RTHT on two dofs
    %
    
    properties
        % input variables
        Hmax
        waiting_time
        time_step
        tol
        max_number_of_iterations
        max_change_rate
        scaling_residual
        noise
        aN;aP
        freq
        % interface variables
        rtc
        channels_tag
        channels_label
        channels_name
        % function
        Residual
        % tensors built here
        LinMat        
        Mn;Kn;Cn
        Mp;Kp;Cp
        Mc;Cc;Kc
        damp_rat
        L2T
        T2N
        S2N
        OrgX
        Zred
        Zred_h1
        Zglobal
        Znt
        Ahglobal
        Rtarget
        Rtarget_h1
        % flags
        live_plots
        record_results
        line_search
        % recorded variables            
        record_residual % while
        record_xsol     % while
        err             % while
        change_rate     % while
        abs_rate     % while
        record_sensors_harm  % once
        record_sensors % once
        invJ            % once
        Xsol            % once#
        Conv
        % figures
        fig40
        fig50
        fig10
    end    
    methods        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%   Main Function   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = BroydenIterations(obj,X0,R0,invJ)
        % Iterations called with: bm = bm.BroydenIterations(X0,R0,invJ)
            it=1;
            obj.err = norm(R0,2);
            DX_0 = -invJ*R0;
            tol0=obj.tol;
            if obj.err(end) <= obj.tol 
                obj.Xsol = X0;
                if obj.record_results == 1
                    obj.record_residual = [obj.record_residual,R0];
                    obj.record_xsol     = [obj.record_xsol,X0];
                end
            end
            while obj.err(end) > obj.tol && it <= obj.max_number_of_iterations
                if obj.rtc.par.stop_all == 0
                    DX_temp = -invJ*R0;
                    obj.change_rate(it) = norm(DX_temp)/norm(DX_0);                    
                    disp('iteration step:')
                    disp(it)
                    disp('change rate:')
                    disp(obj.change_rate(it))                    
                    if obj.change_rate(it) > obj.max_change_rate
                        alpha=obj.change_rate(it)/obj.max_change_rate;
                    else
                        alpha=1;
                    end
                    obj.abs_rate(it) = norm(X0+DX_temp./alpha)/norm(DX_0);
                    %disp('abs rate:')
                    %disp(obj.abs_rate(it))
                    %
                    full_step = - (invJ*R0)./alpha;
                    if obj.line_search && floor((it-1)/10)*10==it-1
                        disp("Line Search")
                        % start line search:                        
                        n_substeps = 3;
                        norm_temp=zeros(n_substeps,1);
                        factor=1/n_substeps:1/n_substeps:1;
                        for iSub = 1:n_substeps
                            Xtemp = X0 + factor(iSub)*full_step;
                            [Rtemp,~] = obj.Residual(obj,Xtemp);
                            norm_temp(iSub)=norm(Rtemp);
                        end
%                         norm_temp_neg=zeros(n_substeps,1);
%                         for iSub = 1:n_substeps
%                             Xtemp = X0 - factor(iSub)*full_step;
%                             [Rtemp,~] = obj.Residual(obj,Xtemp);
%                             norm_temp_neg(iSub)=norm(Rtemp);
%                         end
%                         norm_temp=...
%                             [norm_temp_neg(end:-1:1);norm(R0);norm_temp];
%                         factor=...
%                             [-factor(end:-1:1),0,factor];
                        best = factor(norm_temp==min(norm_temp));
                        disp("best step length")
                        disp(best);
                        worst = factor(norm_temp==max(norm_temp));
                        disp("worst step length")
                        disp(worst);
                        % h=figure(30);%h.WindowState='minimized';
                        % plot(factor,norm_temp)
                        % go back smoothly
                        %for iFactor = 1:-1/n_substeps:best+1/n_substeps
                        %    Xtemp = X0 + iFactor*full_step;
                        %    [~,~] = obj.Residual(obj,Xtemp);
                        %end
                        % end line search
                    else
                        best=1;
                    end
                    obj.Xsol = X0 + best*full_step;
                    [R1,noise_R1] = obj.Residual(obj,obj.Xsol);
                    it=it+1;
                    obj.err(it) = norm(R1);
                    disp('norm of residual:')
                    disp(obj.err(it))
                    disp('std(noise of R):')
                    disp(sqrt(norm(noise_R1)))
                    disp('change in direction of residual:')
                    disp(R1'*R0/(norm(R1)*norm(R0)))
                    disp('tolerance:')
                    disp(tol0)
                    obj.tol = tol0;%min(1.6*sqrt(norm(noise_R1)),tol0);
                    %
                    %
                    obj.fig40.Children.Children(2).XData=1:length(obj.err);
                    obj.fig40.Children.Children(2).YData=obj.err;
                    obj.fig40.Children.Children(1).XData=1:length(obj.err);
                    obj.fig40.Children.Children(1).YData=obj.tol+0*obj.err;

                    if floor(it/20)*20==it && obj.err(end) > 2*obj.tol && obj.err(end)>obj.err(end-1)
                        disp("Numerical Jacobian")
                        e=eye(4);
                        J=e*0;
                        sc=.15;
                        for i=1:4
                            step_dX=max(0.05,sc*((obj.Xsol(i))));
                            [Ri,~]=obj.Residual(obj,obj.Xsol+step_dX*e(:,i));
                            J(:,i) = (Ri-R1)./(step_dX);
                        end
                        
                        invJ =inv(J);
                    else
                        invJ = invJ + ...
                            ((obj.Xsol-X0-invJ*(R1-R0))*(R1-R0)')./norm((R1-R0))^2;
                    end
                    R0 = R1;
                    X0 = obj.Xsol;
                    if obj.record_results == 1
                        obj.record_residual = [obj.record_residual,R0];
                        obj.record_xsol     = [obj.record_xsol,X0]; 
                        obj.invJ     = invJ;
                    end
                else
                    break
                end
            end            
            obj.Conv = obj.err(end) < obj.tol;
            obj.Mn=[];obj.Kn=[];obj.Cn=[];obj.Mp=[];obj.Kp=[];obj.Cp=[];
            % once converged
            if obj.record_results == 1
                N_periods = 2;
                sampling_rate = round(1/obj.time_step);
                obj.rtc.set_stream(0, obj.channels_tag, floor(sampling_rate/obj.freq*N_periods), 0);
                obj.record_sensors=obj.rtc.run_stream(0);
                obj.record_sensors_harm = [...
                    obj.rtc.par.in1_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                    obj.rtc.par.in2_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                    obj.rtc.par.in3_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                    obj.rtc.par.in4_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                    obj.rtc.par.in5_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                    obj.rtc.par.in6_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])'];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%   Common Functions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ImposeVoltage(obj,X)
            obj.rtc.set_par('voltage_freq', obj.freq);
            V7_coef = zeros(16,1);
            V8_coef = zeros(16,1);
            % Impose voltage on shakers and wait
            V7_coef(2:obj.Hmax+1)=X(1:obj.Hmax);
            V7_coef(2+8:obj.Hmax+1+8)=X(obj.Hmax+1:2*obj.Hmax);
            V8_coef(2:obj.Hmax+1)=X(2*obj.Hmax+1:3*obj.Hmax);
            V8_coef(2+8:obj.Hmax+1+8)=X(3*obj.Hmax+1:4*obj.Hmax);            
            % impose it
            obj.rtc.set_par('voltage7_coeffs',V7_coef');
            obj.rtc.set_par('voltage8_coeffs',V8_coef');
            if obj.rtc.par.stop_all ~= 0
                disp("stop!")
                disp(obj.rtc.par.stop_all)
                return
            end
            pause(obj.time_step*obj.waiting_time)
            if obj.live_plots == 1
                LivePlots(obj);
            end
        end
        function X_larger = EnlargeX(obj,H_prev,X_prev)
            % Enlarge the unknowns vector with zeros
            X_larger = [...
                X_prev(1+0*H_prev:1*H_prev);  zeros(obj.Hmax-H_prev,1);...
                X_prev(1+1*H_prev:2*H_prev);  zeros(obj.Hmax-H_prev,1);...
                X_prev(1+2*H_prev:3*H_prev);  zeros(obj.Hmax-H_prev,1);...
                X_prev(1+3*H_prev:4*H_prev);  zeros(obj.Hmax-H_prev,1)...
                ];
        end
        function iJ_larger = EnlargeJ(obj,H_prev,iJ_prev,scA)
            % Enlarge the jacobian matrix
            iJ_larger = eye(4*obj.Hmax)/scA;
            ind_prev = [0*obj.Hmax+(1:H_prev),...
                        1*obj.Hmax+(1:H_prev),...
                        2*obj.Hmax+(1:H_prev),...
                        3*obj.Hmax+(1:H_prev)];
            iJ_larger(ind_prev,ind_prev) = iJ_prev;
        end
        function X_comp = ComposeX(obj,H_prev,X_prev,X_higher)
            % Enlarge the unknowns vector with a higher harmonic
            X_comp = [...
                X_prev(1+0*H_prev:1*H_prev);  X_higher(1+0*1:1*1);...
                X_prev(1+1*H_prev:2*H_prev);  X_higher(1+1*1:2*1);...
                X_prev(1+2*H_prev:3*H_prev);  X_higher(1+2*1:3*1);...
                X_prev(1+3*H_prev:4*H_prev);  X_higher(1+3*1:4*1)...
                ];
        end
        function R_higher = SelectR(obj,H_prev,R_prev)
            % Selects the highest harmonics of the residual
            R_higher = [...
                R_prev(1*H_prev);...
                R_prev(2*H_prev);...
                R_prev(3*H_prev);...
                R_prev(4*H_prev)...
                ];
        end
        function obj = AssembleRtarget(obj,FirstHarmTarget)
            % Target displacement inner loop style
            obj.Rtarget=[...
                FirstHarmTarget(1);zeros(obj.Hmax-1,1);...
                FirstHarmTarget(2);zeros(obj.Hmax-1,1);...
                FirstHarmTarget(3);zeros(obj.Hmax-1,1);...
                FirstHarmTarget(4);zeros(obj.Hmax-1,1)];
        end
        function LivePlots(obj)
            N_periods = 2;
            sampling_rate = round(1/obj.time_step);
            periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            obj.rtc.set_stream(0, obj.channels_tag, floor(sampling_rate/obj.freq*N_periods), 0);
            st=obj.rtc.run_stream(0);
            
            N_periods = N_periods-1;
            Fourier_periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            
            iDFT = [sin(2*pi*(1:obj.Hmax).*Fourier_periods'),...
                cos(2*pi*(1:obj.Hmax).*Fourier_periods')];
            
            coeff = [...
                obj.rtc.par.in1_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in2_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in3_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in4_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in5_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in6_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])'];
            
            if obj.fig10==0
                obj.fig10=figure(10);
                t=tiledlayout(4,2);
                for iCh=1:6
                    ax(iCh) = nexttile;
                    plot(Fourier_periods,...
                        iDFT*coeff(:,iCh),'r.'...
                        );
                    hold on;
                    const = (max(st(iCh,:))+min(st(iCh,:)))/2;
                    plot(st(9,1)/(2*pi)+periods-1,(st(iCh,:)-const),'b');
                    hold off
                    xlim([0 1])
                    ylabel(obj.channels_label{iCh});
                end
                for iCh=7:8
                    ax(iCh) = nexttile;
                    const = (max(st(iCh,:))+min(st(iCh,:)))/2;
                    plot(st(9,1)/(2*pi)+periods-1,(st(iCh,:)-const),'b');
                    xlim([0 1])
                    ylabel(obj.channels_label{iCh});
                end          
                linkaxes(ax,'x');
            else
                for iCh=1:6
                    obj.fig10.Children.Children(9-iCh).Children(2).XData=Fourier_periods;
                    obj.fig10.Children.Children(9-iCh).Children(2).YData=iDFT*coeff(:,iCh);
                    const = (max(st(iCh,:))+min(st(iCh,:)))/2;
                    obj.fig10.Children.Children(9-iCh).Children(1).XData=st(9,1)/(2*pi)+periods-1;
                    obj.fig10.Children.Children(9-iCh).Children(1).YData=st(iCh,:)-const;
                end
                for iCh=7:8
                    obj.fig10.Children.Children(9-iCh).Children(1).XData=st(9,1)/(2*pi)+periods-1;
                    obj.fig10.Children.Children(9-iCh).Children(1).YData=st(iCh,:);
                end
            end
        end
        function NonLivePlots(obj)
            N_periods = 2;
            sampling_rate = round(1/obj.time_step);
            periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            obj.rtc.set_stream(0, obj.channels_tag, floor(sampling_rate/obj.freq*N_periods), 0);
            st=obj.rtc.run_stream(0);
            
            N_periods = N_periods-1;
            Fourier_periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            
            iDFT = [sin(2*pi*(1:obj.Hmax).*Fourier_periods'),...
                cos(2*pi*(1:obj.Hmax).*Fourier_periods')];
            
            coeff = [...
                obj.rtc.par.in1_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in2_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in3_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in4_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in5_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in6_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])'];
            % figure(1)
            title('v   Front  |    Back    v')
            t = tiledlayout(4,2);
            for iCh=1:6
                ax(iCh) = nexttile;
                plot(Fourier_periods,...
                    iDFT*coeff(:,iCh),'r'...
                    ,'DisplayName',strcat(obj.channels_name{iCh}," ",'Fourier')...
                    );
                hold on;
                const = (max(st(iCh,:))+min(st(iCh,:)))/2;
                plot(st(9,1)/(2*pi)+periods-1,(st(iCh,:)-const),'b',...
                'DisplayName',strcat(obj.channels_name{iCh}," ",'Real signal'));
                xlabel('periods')
                ylabel(obj.channels_label{iCh});
                legend();
            end
            for iCh=7:8
                ax(iCh) = nexttile;
                const = (max(st(iCh,:))+min(st(iCh,:)))/2;
                plot(st(9,1)/(2*pi)+periods-1,(st(iCh,:)-const),'b',...
                'DisplayName',strcat(obj.channels_name{iCh}," ",'Real signal'));
                xlabel('periods')
                ylabel(obj.channels_label{iCh});
                legend();
            end            
            linkaxes(ax,'x');
        end
        function OfflinePlotsFourier(obj,X,name)
            N_periods = 1;
            sampling_rate = round(1/obj.time_step);
            Fourier_periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            iDFT = [sin(2*pi*(1:obj.Hmax).*Fourier_periods'),...
                    cos(2*pi*(1:obj.Hmax).*Fourier_periods')];
            if name(1)=='V'
                plot(Fourier_periods,iDFT*X,'-','LineWidth',2,'DisplayName',name,'Color',[0,.7,.4])
            elseif name(1)=='P'
                plot(Fourier_periods,iDFT*X,'--','LineWidth',2,'DisplayName',name,'Color',[0,0.4,.8])
            else
                plot(Fourier_periods,iDFT*X,'LineWidth',2,'DisplayName',name)
            end
            
            hold on;
        end
        function OfflinePlotsTime(obj,X,lag,name)
            N_periods=2;
            sampling_rate = round(1/obj.time_step);
            periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            const = (max(X)+min(X))/2;
            plot(lag/(2*pi)+periods-1,(X-const),'k','LineWidth',.5,'DisplayName',name)
            xlim([0,periods(floor(end/2))])
            hold on;
        end
        function [st,coeff] = RunStream(obj,N_periods)
            sampling_rate = round(1/obj.time_step);
            %floor(sampling_rate/obj.freq*N_periods)
            obj.rtc.set_stream(0, obj.channels_tag, floor(sampling_rate/obj.freq*N_periods), 0);
            st=obj.rtc.run_stream(0);
            coeff = [...
                obj.rtc.par.in1_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in2_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in3_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in4_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in5_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])',...
                obj.rtc.par.in6_coeffs([2:obj.Hmax+1,10:obj.Hmax+9])'];
        end
        function PlotStream(obj,N_periods,st,coeff)
            sampling_rate = round(1/obj.time_step);
            periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            % obj.rtc.set_stream(0, obj.channels_tag, floor(sampling_rate/obj.freq*N_periods), 0);
            % st=obj.rtc.run_stream(0);
            
            N_periods = N_periods-1;
            Fourier_periods = [1:floor(sampling_rate/obj.freq*N_periods)]/floor(sampling_rate/obj.freq);
            
            iDFT = [sin(2*pi*(1:obj.Hmax).*Fourier_periods'),...
                cos(2*pi*(1:obj.Hmax).*Fourier_periods')];
            
            
            figure(1)
            %title('v   Front  |    Back    v')
            t = tiledlayout(4,2);
            for iCh=1:6
                ax(iCh) = nexttile;
                plot(1 + Fourier_periods,...
                    iDFT*coeff(:,iCh),'r'...
                    );
                    %,'DisplayName',strcat(obj.channels_name{iCh}," ",'Fourier'));
                hold on;
                const = (max(st(iCh,:))+min(st(iCh,:)))/2;
                plot(st(9,1)/(2*pi)+periods,(st(iCh,:)-const),'b');
                %'DisplayName',strcat(obj.channels_name{iCh}," ",'Real signal'));
                %xlabel('periods')
                ylabel(obj.channels_label{iCh});
                %legend();
            end
            for iCh=7:8
                ax(iCh) = nexttile;
                const = (max(st(iCh,:))+min(st(iCh,:)))/2;
                plot(st(9,1)/(2*pi)+periods,(st(iCh,:)-const),'b');
                %'DisplayName',strcat(obj.channels_name{iCh}," ",'Real signal'));
                %xlabel('periods')
                ylabel(obj.channels_label{iCh});
                %legend();
            end            
            linkaxes(ax,'x');
        end
        function obj = InterfaceInit(obj,lab)
            if lab
                obj.rtc=it_2out_num_interface();
                obj.rtc.set_par('voltage7_coeffs',zeros(1,16));
                obj.rtc.set_par('voltage8_coeffs',zeros(1,16));
                obj.rtc.par.stop_all=0;
                obj.rtc.set_par('voltage_freq', obj.freq);
            else
                obj.rtc={};
                obj.rtc.par={};
                obj.rtc.par.stop_all = 0;
            end
            
            input_channels_tag = {'in1','in2','in3','in4','in5','in6'};
            output_channels_tag = {'voltage7','voltage8'};
            time_channels_tag = {'time_mod_2pi'};
            input_channels_name = {'acc front','acc back',...
                'force front','force back',...% Range: 2.224 kN, Sensitivity: 2248 mV/kN, Gain: 100
                'laser front','laser back'};% Range 10 mm, Sensitivity: high, Gain:
            output_channels_name = {'shaker front','shaker back'};
            
            obj.channels_label = {'mm/s^2','mm/s^2','N','N','mm','mm','V','V'};
            obj.channels_tag = [input_channels_tag,output_channels_tag,time_channels_tag];
            obj.channels_name = [input_channels_name,output_channels_name];
            
            obj.fig50=figure(50);plot(0:1);hold on;plot(1:2,'--');hold off;
            obj.fig40=figure(40);semilogy(1:2);hold on;semilogy(1:2,'.-');hold off;
            obj.fig10=figure(10);t=tiledlayout(4,2);
            for iCh=1:6;ax(iCh) = nexttile;plot(0:0.5:1,0:0.5:1,'r.');hold on;plot(0:1,0:1,'b');hold off;xlim([0 1]);ylabel(obj.channels_label{iCh});end
            for iCh=7:8;ax(iCh) = nexttile;plot(0:1,0:1,'b');xlim([0 1]);ylabel(obj.channels_label{iCh});end;linkaxes(ax,'x');

        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%   Residual outer loop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = ComputeDynamicMatricesNumerical(obj)
            % provides KN, CN, MN with modal damping
            n = 500; % number of elements is seldom changed
            % best number to scale mom and rot is aN:
            dl = 1/n; % such dl makes: rot = rot*aN, mom = mom/aN
            % [Force]=N, [Time]=s, [Length]=mm
            % elemental matrices non dimensional
            Ke = [...
                12,     6*dl,       -12,          6*dl;...
                6*dl,   4*dl^2,     -6*dl,        2*dl^2;...
                -12,    -6*dl,        12,         -6*dl;...
                6*dl,   2*dl^2,     -6*dl,        4*dl^2];            
            Me = 1/420*[...
                156,        22*dl,      54,     -13*dl;...
                22*dl,      4*dl^2,     13*dl,  -3*dl^2;...
                54,         13*dl,      156,    -22*dl;...
                -13*dl,    -3*dl^2,    -22*dl,   4*dl^2];            
            KN = zeros(2*n);
            MN = zeros(2*n);
            for iDof = 1:n-1
                KN( 1+2*(iDof-1) : 2*(iDof+1), 1+2*(iDof-1): 2*(iDof+1) ) = ...
                    KN( 1+2*(iDof-1) : 2*(iDof+1), 1+2*(iDof-1) : 2*(iDof+1) ) + ...
                    Ke;
            end
            for iDof = 1:n-1
                MN( 1+2*(iDof-1) : 2*(iDof+1), 1+2*(iDof-1): 2*(iDof+1) ) = ...
                    MN( 1+2*(iDof-1) : 2*(iDof+1), 1+2*(iDof-1) : 2*(iDof+1) ) + ...
                    Me;
            end
            % parameters of numerical model
            rho=8.03E-9*(6.05/6.4789)^-2;%
            E=217E3;% 216.3E3 N/mm^2
            I = 25.4*(1)^3/12 ; % mm^4
            A = 25.4*(1); % mm^2
            % physical matrices (now the aP is used)
            obj.Kp = (n/obj.aP)^3*E*I*KN;
            obj.Mp = (obj.aP/n)*rho*A*MN;
            % clamped root
            KN = KN(3 :end, 3 :end);
            MN = MN(3 :end, 3 :end);
            % numerical matrices (now the aN is used)
            obj.Kn = (n/obj.aN)^3*E*I*KN;
            obj.Mn = (obj.aN/n)*rho*A*MN;
            % build modal damping matrix
            [V,D]=eigs(obj.Kn,obj.Mn,2*n-2,'smallestabs');
            Xi = 2*obj.damp_rat*sqrt(D);
            invV = inv(V);
            obj.Cn = invV'*Xi*invV;
            [V,D]=eigs(obj.Kp,obj.Mp,2*n,'smallestabs');
            Xi = 2*obj.damp_rat*sqrt(D);
            invV = inv(V);
            obj.Cp = invV'*Xi*invV;
            % check eigs:
            % sqrt(eigs(obj.Kn,obj.Mn,3,'smallestabs'))'/2/pi
            % ([1.875,4.694091132974175,7.855]).^2./(2*pi*aN^2).*sqrt(E*I/rho/A)
            % check static displ:
            % F_tip = 1;
            % F = zeros(2*(n-1),1);
            % F(end-1) = F_tip;
            % X = inv(obj.Kn)*F;
            % X(end-1)
            % F_tip*aN^3/(3*E*I)
            %
            %
            % identified clamp matrices in right dofs organisation:
            load('identified_matrices_HighFreqNew.mat',...
                               'k11','k12','k21','k22',...
                               'c11','c12','c21','c22',...
                               'm11','m12','m21','m22');
            % Boundary dofs reorganised
            % Znt*[disS; disC; rotS; rotC] = 
            %          [forS; forC; momS; momC]
            obj.Kc = [  k11, 0, k12, 0;...
                        0, k11, 0, k12;...
                        k21, 0, k22, 0;...
                        0, k21, 0, k22];
            %
            obj.Cc = [  0,-c11, 0,-c12;...
                        c11, 0, c12, 0;...                        
                        0,-c21, 0,-c22
                        c21, 0, c22, 0;];
            %
            obj.Mc =  [  m11, 0, m12, 0;...
                        0, m11, 0, m12;...
                        m21, 0, m22, 0;...
                        0, m21, 0, m22];
        end 
        function obj = ComputeRotationMatrices(obj)
            % Reduced on boundary dofs reorganised
            % Znt*[disS; disC; rotS; rotC] = 
            %          [forS; forC; momS; momC]
            %
            % [disS; rotS*aN; disC; rotC*aN] = OrgX*[disS; disC; rotS; rotC]
            obj.OrgX = [1,0,0,0; 0,0,obj.aN,0; 0,1,0,0; 0,0,0,obj.aN];
            %
            % [forS; forC; momS; momC] = OrgF*[forS; momS/aN; forC; momC/aN]
            % OrgF = [1,0,0,0; 0,0,1,0; 0,obj.aN,0,0; 0,0,0,obj.aN];
            % OrgF = OrgX';
            % 
            % OrgF*Zn*OrgX*[disS; disC; rotS; rotC] = 
            %          [forS; forC; momS; momC]
            %
            % from numerical model to sensors:
            L_length = 40;
            T_length = 111.6;
            % from lasers to transducers:
            % #obj.L2T = [eye(2),                  zeros(2);...
            % #      T_length/L_length*eye(2),(T_length-L_length)/L_length*eye(2)];
            obj.L2T = [eye(2),                  zeros(2);...
                   -(T_length-L_length)/L_length*eye(2),T_length/L_length*eye(2),];
            % from transducers to numerical:
            % obj.T2N'*[forS; forC; momS; momC] = 
            %         [f1S; f1C; f2S; f2C];
            % obj.T2N*[u1S; u1C; u2S; u2C] = 
            %         [disS; disC; momS; momC];
            obj.T2N = [eye(2),                  zeros(2);...
                   eye(2)/T_length,        -eye(2)/T_length];
            %
            % from strain gauges to numerical:
            % obj.S2N*[sg1S; sg1C; sg2S; sg2C] = 
            %         [forS; forC; momS; momC];
            c1=-437.1645;%437.1645*2.21;3074.15 % calibration of SG to M
            c2=-462.8841;%462.8841*2.185;3640.36 % calibration of SG to M
            d_rel=6.557; % distance between SG1 and SG2
            d0=7; % distance between SG1 and clamp
            % all in the box now!
            % th1=-0.035;
            th=-0.03;
            % R=[cos(th1),sin(th1);-sin(th1),cos(th1)];
            R=[cos(th),sin(th);-sin(th),cos(th)];
            obj.S2N = [ R,   0*R;...
                        0*R,   R];
%             obj.S2N = [ c1/d_rel*eye(2),   -c2/d_rel*eye(2);...
%                          c1*eye(2),  0*eye(2)];
            %
        end
        function obj = ComputeReducedDynamic(obj,om)     
            % provide Zred: reduced dynamical stiffness matrix at boundary
            % Zn*[disS; rotS*aN; disC; rotC*aN] = 
            %          [forS; momS/aN; forC; momC/aN]
            % Kn,Mn,Cn are defined once and for all
            DN = -om^2*obj.Mn + obj.Kn;
            ON = om*obj.Cn;
            % Dynamic stiffness of the whole numerical structure for one harmonic
            % ZN =  [  DN  ,   -ON ;...
            %          ON  ,    DN ];
            % Dofs organisation
            % ZN*[Xs;Xc] = [Fs;Fc];  
            % X = [dis aN*rot dis aN*rot ...]
            % F = [for mom/aN for mom/aN ...]
            % Internal dofs and boundary dofs matrices
            Zii =  [  DN(1:end-2,  1:end-2)  , -ON(1:end-2,  1:end-2)  ; ...
                     +ON(1:end-2,  1:end-2)  ,  DN(1:end-2,  1:end-2)  ];
            Zib =  [  DN(1:end-2,  end-1:end), -ON(1:end-2,  end-1:end); ...
                     +ON(1:end-2,  end-1:end),  DN(1:end-2,  end-1:end)];
            Zbi =  [  DN(1:end-2,  end-1:end), +ON(1:end-2,  end-1:end); ...
                     -ON(1:end-2,  end-1:end),  DN(1:end-2,  end-1:end)]';
            Zbb =  [  DN(end-1:end,end-1:end), -ON(end-1:end,end-1:end); ...
                     +ON(end-1:end,end-1:end),  DN(end-1:end,end-1:end)];
            % Reduced on boundary dofs
            % Zn*[disS; rotS*aN; disC; rotC*aN] = 
            %          [forS; momS/aN; forC; momC/aN]
            Zn = Zbb - Zbi*(Zii\Zib);
            % Reduced on boundary dofs reorganised
            % Znt*[disS; disC; rotS; rotC] = 
            %          [forS; forC; momS; momC]
            obj.Znt = obj.OrgX'*Zn*obj.OrgX;
            % 
            % with identified clamp
            obj.Zred = (obj.T2N'*obj.Znt*obj.T2N)*obj.L2T +...
                -(-om^2*obj.Mc +om*obj.Cc +obj.Kc)*obj.L2T;
            
            % without:
            % obj.Zred = (obj.T2N'*Znt*obj.T2N)*obj.L2T;
        end
        function obj = AssembleAhglobal(obj,Ex_h1)
            obj = ComputeDynamicMatricesNumerical(obj);
            obj = ComputeRotationMatrices(obj);
            obj.Zglobal = zeros(4*obj.Hmax);            
            for h = 1:obj.Hmax
                % here Zred is build for each harmonic
                obj = ComputeReducedDynamic(obj,h*(obj.freq*2*pi));
                % Zred is the single harm dynamic stiffness matrix 
                % that relates laser input to forces at transducers location
                % it is fed into Zglobal
                ind = h+[0,1,2,3]*obj.Hmax;
                obj.Zglobal(ind,ind) = obj.Zred;
                % Excitation vector on first harm, sinusoidal force at interface:
                % Ex_h1 = [.1;0;0;0];
                if h==1
                    obj.Zred_h1=obj.Zred;
                    obj.Rtarget_h1 = -obj.T2N'*Ex_h1;
                end
            end
            obj = AssembleRtarget(obj,obj.Rtarget_h1);
            % forces on clamp from single computation
            % load('identified_forces.mat','F');
            % obj.Rtarget=obj.Rtarget + ...
            %     [F([0+1:0+obj.Hmax,length(F)/2+1:length(F)/2+obj.Hmax],1);...
            %    F([0+1:0+obj.Hmax,length(F)/2+1:length(F)/2+obj.Hmax],2)];
        end
        function [R,noise_R] = BuildResidualOuterLoop(obj)
            % Residual outer loop style
            % [TransFS; TransFC; TransBS; TransBC] + Rtarget 
            %        + Zglobal*[LaserFS; LaserBS; LaserFC; LaserBC]
            %           ;                      
            R=obj.Zglobal*[...
                obj.rtc.par.in5_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_ave(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(8+2:8+obj.Hmax+1)';...
                ] + ...
                [...
                obj.rtc.par.in3_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in3_coeffs_ave(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_ave(8+2:8+obj.Hmax+1)';...
                ] + ...
                obj.Rtarget;% includes forces to clamp from single computation
            noise_R=obj.Zglobal*[...
                obj.rtc.par.in5_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_var(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_var(0+2:0+obj.Hmax+1)';...                
                obj.rtc.par.in6_coeffs_var(8+2:8+obj.Hmax+1)';...
                ] + ...
                [...
                obj.rtc.par.in3_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in3_coeffs_var(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_var(8+2:8+obj.Hmax+1)';...
                ];
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
                noise_R=noise_R.*obj.scaling_residual;
            end
                        
        end
        function [R,noise_R] = ResidualOuterLoop(obj,X)
            n_periods=10;
            R=zeros(size(X,1),n_periods);
            for i=1:n_periods
                ImposeVoltage(obj,X);
                [R(:,i),noise_R] = BuildResidualOuterLoop(obj);
            end
            %
            std=norm(R(:,end-1)-R(:,end));
            while std(end) > sqrt(norm(noise_R))/2
                n_periods=n_periods+1;
                ImposeVoltage(obj,X);
                [R_new,noise_R] = BuildResidualOuterLoop(obj);
                R=[R,R_new];
                std=[std,norm(R(:,end-1)-R(:,end))];
                if obj.fig50==0
                    obj.fig50=figure(50);
                end
                obj.fig50;
                plot(std);hold on
                plot(std*0+sqrt(norm(noise_R))/2);hold off
            end
            R=R(:,end);
            %disp("transient decayed in cycles")
            %disp(n_periods)
            
        end
        function R   = ResidualOuterLoop2(obj,X)
            % Residual outer loop style only on first harm
            % Remaining res is just setting displ higher harm to zero
            % [TransFS; TransFC; TransBS; TransBC] + Rtarget 
            %        - Zglobal*[LaserFS; LaserBS; LaserFC; LaserBC]
            %           ;            
            ImposeVoltage(obj,X);           
            R=[...
                obj.rtc.par.in5_coeffs_ave(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in5_coeffs_ave(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in6_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(8+2:8+obj.Hmax+1)'...
                ];
            noise_R = [...
                obj.rtc.par.in5_coeffs_var(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in5_coeffs_var(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in6_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_var(8+2:8+obj.Hmax+1)'...
                ];
            ind_h1=[1,obj.Hmax+1,2*obj.Hmax+1,3*obj.Hmax+1];
            R(ind_h1)=-obj.Zred_h1*[...
                obj.rtc.par.in5_coeffs_ave( 2)';...
                obj.rtc.par.in5_coeffs_ave(10)';...
                obj.rtc.par.in6_coeffs_ave( 2)';...
                obj.rtc.par.in6_coeffs_ave(10)';...
                ] + ...
                [...
                obj.rtc.par.in3_coeffs_ave( 2)';...
                obj.rtc.par.in3_coeffs_ave(10)';...
                obj.rtc.par.in4_coeffs_ave( 2)';...
                obj.rtc.par.in4_coeffs_ave(10)';...
                ] + ...
                obj.Rtarget_h1;
            noise_R(ind_h1)=-obj.Zred_h1*[...
                obj.rtc.par.in5_coeffs_var( 2)';...
                obj.rtc.par.in5_coeffs_var(10)';...
                obj.rtc.par.in6_coeffs_var( 2)';...
                obj.rtc.par.in6_coeffs_var(10)';...
                ] + ...
                [...
                obj.rtc.par.in3_coeffs_var( 2)';...
                obj.rtc.par.in3_coeffs_var(10)';...
                obj.rtc.par.in4_coeffs_var( 2)';...
                obj.rtc.par.in4_coeffs_var(10)';...
                ];
            disp('std(noise of R):')
            disp(sqrt(norm(noise_R)))
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%   Residual outer loop  SG  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = ComputeReducedDynamicSG(obj,om)     
            % provide Zred: reduced dynamical stiffness matrix at boundary
            % Zn*[disS; rotS*aN; disC; rotC*aN] = 
            %          [forS; momS/aN; forC; momC/aN]
            % Kn,Mn,Cn are defined once and for all
            DN = -om^2*obj.Mn + obj.Kn;
            ON = om*obj.Cn;
            % Dynamic stiffness of the whole numerical structure for one harmonic
            % ZN =  [  DN  ,   -ON ;...
            %          ON  ,    DN ];
            % Dofs organisation
            % ZN*[Xs;Xc] = [Fs;Fc];  
            % X = [dis aN*rot dis aN*rot ...]
            % F = [for mom/aN for mom/aN ...]
            % Internal dofs and boundary dofs matrices
            Zii =  [  DN(1:end-2,  1:end-2)  , -ON(1:end-2,  1:end-2)  ; ...
                     +ON(1:end-2,  1:end-2)  ,  DN(1:end-2,  1:end-2)  ];
            Zib =  [  DN(1:end-2,  end-1:end), -ON(1:end-2,  end-1:end); ...
                     +ON(1:end-2,  end-1:end),  DN(1:end-2,  end-1:end)];
            Zbi =  [  DN(1:end-2,  end-1:end), +ON(1:end-2,  end-1:end); ...
                     -ON(1:end-2,  end-1:end),  DN(1:end-2,  end-1:end)]';
            Zbb =  [  DN(end-1:end,end-1:end), -ON(end-1:end,end-1:end); ...
                     +ON(end-1:end,end-1:end),  DN(end-1:end,end-1:end)];
            % Reduced on boundary dofs
            % Zn*[disS; rotS*aN; disC; rotC*aN] = 
            %          [forS; momS/aN; forC; momC/aN]
            Zn = Zbb - Zbi*(Zii\Zib);
            % Reduced on boundary dofs reorganised
            % Znt*[disS; disC; rotS; rotC] = 
            %          [forS; forC; momS; momC]
            obj.Znt = obj.OrgX'*Zn*obj.OrgX;
            % 
            % 
            obj.Zred = (obj.T2N'*obj.Znt*obj.T2N)*obj.L2T;
            
        end
        function obj = AssembleAhglobalSG(obj,Ex_h1)
            obj = ComputeDynamicMatricesNumerical(obj);
            obj = ComputeRotationMatrices(obj);
            obj.Zglobal = zeros(4*obj.Hmax);            
            for h = 1:obj.Hmax
                % here Zred is build for each harmonic
                obj = ComputeReducedDynamicSG(obj,h*(obj.freq*2*pi));
                % Zred is the single harm dynamic stiffness matrix 
                % that relates laser input to forces at transducers location
                % it is fed into Zglobal
                ind = h+[0,1,2,3]*obj.Hmax;
                obj.Zglobal(ind,ind) = obj.Zred;
                % Excitation vector on first harm, sinusoidal force at interface:
                % Ex_h1 = [.1;0;0;0];
                if h==1
                    obj.Zred_h1=obj.Zred;
                    obj.Rtarget_h1 = -obj.T2N'*Ex_h1;
                end
            end
            obj = AssembleRtarget(obj,obj.Rtarget_h1);
            % forces on clamp from single computation
            % load('identified_forces.mat','F');
            % obj.Rtarget=obj.Rtarget + ...
            %     [F([0+1:0+obj.Hmax,length(F)/2+1:length(F)/2+obj.Hmax],1);...
            %    F([0+1:0+obj.Hmax,length(F)/2+1:length(F)/2+obj.Hmax],2)];
        end
        function [R,noise_R,D,noise_D] = BuildResidualSG(obj)
            % Residual outer loop style
            % [TransFS; TransFC; TransBS; TransBC] + Rtarget 
            %        + Zglobal*[LaserFS; LaserBS; LaserFC; LaserBC]
            %           ;      
            D=[...
                obj.rtc.par.in5_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_ave(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(8+2:8+obj.Hmax+1)'];
            noise_D=[...
                obj.rtc.par.in5_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_var(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_var(0+2:0+obj.Hmax+1)';...                
                obj.rtc.par.in6_coeffs_var(8+2:8+obj.Hmax+1)'];
            
            R=obj.Zglobal*[...
                obj.rtc.par.in5_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_ave(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(8+2:8+obj.Hmax+1)';...
                ] + ...
                obj.S2N*[...%obj.T2N'* in the box now!
                obj.rtc.par.in1_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in1_coeffs_ave(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in2_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in2_coeffs_ave(8+2:8+obj.Hmax+1)';...
                ] + ...
                obj.Rtarget;% includes forces to clamp from single computation
            noise_R=obj.Zglobal*[...
                obj.rtc.par.in5_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_var(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_var(0+2:0+obj.Hmax+1)';...                
                obj.rtc.par.in6_coeffs_var(8+2:8+obj.Hmax+1)';...
                ] + ...
                obj.S2N*[...
                obj.rtc.par.in1_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in1_coeffs_var(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in2_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in2_coeffs_var(8+2:8+obj.Hmax+1)';...
                ];
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
                noise_R=noise_R.*obj.scaling_residual;
            end
                        
            D=R;noise_D=noise_R;
        end
        function [R,noise_R] = ResidualSG(obj,X)
            n_periods=30;
            D=zeros(size(X,1),n_periods);
            for i=1:n_periods
                ImposeVoltage(obj,X);
                [R,noise_R,D(:,i),noise_D] = BuildResidualSG(obj);
            end
            %
            std=norm(D(:,end-1)-D(:,end));
            max_std=.001;
            while std(end) > max_std%.15*sqrt(norm(noise_D))
                n_periods=n_periods+1;
                ImposeVoltage(obj,X);
                [R,noise_R,D_new,~] = BuildResidualSG(obj);
                D=[D,D_new];
                std=[std,norm(D(:,end-1)-D(:,end))];
                %
                %
                obj.fig50.Children.Children(2).XData=1:length(std);
                obj.fig50.Children.Children(2).YData=std;
                obj.fig50.Children.Children(1).XData=1:length(std);
                obj.fig50.Children.Children(1).YData=std*0+max_std;
            end
            %R=R(:,end);
            %disp("transient decayed in cycles")
            %disp(n_periods)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%   Residual inner loop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [R,noise_R] = ResidualInnerLoopDispl(obj,X)
            % Residual inner loop style
            ImposeVoltage(obj,X);
            R=obj.Rtarget-[...
                obj.rtc.par.in5_coeffs_ave(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in5_coeffs_ave(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in6_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(8+2:8+obj.Hmax+1)'...
                ];
            noise_R = -[...
                obj.rtc.par.in5_coeffs_var(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in5_coeffs_var(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in6_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_var(8+2:8+obj.Hmax+1)'...
                ];
            disp('std(noise of R):')
            disp(sqrt(norm(noise_R)))
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
            end
        end
        function [R,noise_R] = ResidualInnerLoopForce(obj,X)
            % Residual inner loop style
            ImposeVoltage(obj,X);
            R=obj.Rtarget-[...
                obj.rtc.par.in3_coeffs_ave(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in3_coeffs_ave(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in4_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_ave(8+2:8+obj.Hmax+1)'...
                ];
            noise_R = -[...
                obj.rtc.par.in3_coeffs_var(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in3_coeffs_var(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in4_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_var(8+2:8+obj.Hmax+1)'...
                ];
            disp('std(noise of R):')
            disp(sqrt(norm(noise_R)))
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
            end
        end
        function [R,noise_R] = ResidualInnerLoopAcc(obj,X)
            % Residual inner loop style
            ImposeVoltage(obj,X);
            R=obj.Rtarget-[...
                obj.rtc.par.in1_coeffs_ave(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in1_coeffs_ave(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in2_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in2_coeffs_ave(8+2:8+obj.Hmax+1)'...
                ];
            noise_R = -[...
                obj.rtc.par.in1_coeffs_var(0+2:0+obj.Hmax+1)';...sin
                obj.rtc.par.in1_coeffs_var(8+2:8+obj.Hmax+1)';...cos
                obj.rtc.par.in2_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in2_coeffs_var(8+2:8+obj.Hmax+1)'...
                ];
            disp('std(noise of R):')
            disp(sqrt(norm(noise_R)))
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%   Residual to numerically test   %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [R,noise_R] = ResidualTest(obj,X)
            load('LinMatTest.mat','A','B');
            obj.LinMat=A(1:4*obj.Hmax,1:4*obj.Hmax);
            obj.Rtarget=B(1:4*obj.Hmax);
            % simulate multiple call of residual with noise
            n_periods=3;
            R=obj.Rtarget - obj.LinMat*X + rand(1)*obj.noise;
            
            i=1;
            max_std=.002;std=max_std*1.1;
            while std(end) > max_std || i < n_periods%.15*sqrt(norm(noise_D))
                    i=i+1;
                    R_new = obj.Rtarget - obj.LinMat*X + rand(1)*obj.noise + obj.noise*30/i;
                    R=[R,R_new];
                    std=[std,norm(R(:,end-1)-R(:,end))];
                    %                   
                    obj.fig50.Children.Children(2).XData=1:length(std);
                    obj.fig50.Children.Children(2).YData=std;
                    obj.fig50.Children.Children(1).XData=1:length(std);
                    obj.fig50.Children.Children(1).YData=std*0+max_std;
                    pause(.05)
            end
            
            R=R(:,end);
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
            end
            noise_R=obj.noise;
            pause(.15)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%   Old functions for residual   %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = ComputeReducedDynamicOld(obj,om)     
            % provide Zred: reduced dynamical stiffness matrix at boundary
            % Zred*[disS; rotS*aN; disC; rotC*aN] = 
            %          [forS; momS/aN; forC; momC/aN]
            DN = -om^2*obj.Mn + obj.Kn;
            ON = om*obj.Cn;
            % Dynamic stiffness for one harmonic
            % ZN =  [  DN  ,   ON ;...
            %         -ON  ,   DN ];
            % Dofs organisation
            % ZN*[Xs;Xc] = [Fs;Fc];  
            % X = [dis aN*rot dis aN*rot ...]
            % F = [for mom/aN for mom/aN ...]
            % Internal dofs and boundary dofs matrices
            Zii =  [  DN(1:end-2,  1:end-2)  ,  ON(1:end-2,  1:end-2)  ; ...
                     -ON(1:end-2,  1:end-2)  ,  DN(1:end-2,  1:end-2)  ];
            Zib =  [  DN(1:end-2,  end-1:end),  ON(1:end-2,  end-1:end); ...
                     -ON(1:end-2,  end-1:end),  DN(1:end-2,  end-1:end)];
            Zbb =  [  DN(end-1:end,end-1:end),  ON(end-1:end,end-1:end); ...
                     -ON(end-1:end,end-1:end),  DN(end-1:end,end-1:end)];
            % Reduced on boundary dofs
            % Zred*[disS; rotS*aN; disC; rotC*aN] = 
            %          [forS; momS/aN; forC; momC/aN]
            obj.Zred = Zbb - Zib'*(Zii\Zib);            
        end
        function obj = AssembleAhglobalOld(obj,Ex_h1)
            % From beam dof to measured dofs
            % for each harmonic:
            % [disS; rotS*aN; disC; rotC*aN] = inv(Zred)*[forS; momS/aN; forC; momC/aN]
            % [LaserFS; LaserBS; LaserFC; LaserBC] = Ld*[disS; rotS*aN; disC; rotC*aN];
            % [forS; momS/aN; forC; momC/aN] = fT*[TransFS; TransBS; TransFC; TransBC];
            % [LaserFS; LaserBS; LaserFC; LaserBC] = Ld*inv(Zred)*fT*[TransFS; TransBS; TransFC; TransBC];
            % A_h = Ld*inv(Zred)*fT;
            % [LaserFS; LaserBS; LaserFC; LaserBC] = A_h*[TransFS; TransBS; TransFC; TransBC];
            L_length = 40;
            T_length = 111.6;
            Ld = [1,0,0,0;    1,-L_length/obj.aN,0,0;  0,0,1,0;    0,0,1,-L_length/obj.aN];
            fT = [-1,-1,0,0;  0,+obj.aN/T_length,0,0;  0,0,-1,-1;  0,0,0,+obj.aN/T_length];
            % Excitation vector on first harm, sinusoidal force at interface:
            % Ex_h1 = [.1;0;0;0];
            % here Kn,Mn,Cn are defined once and for all
            obj = ComputeDynamicMatricesNumerical(obj);
            obj.Ahglobal = zeros(4*obj.Hmax);            
            for h = 1:obj.Hmax
                % here Zred is build for each harmonic
                obj = ComputeReducedDynamic(obj,h*(obj.freq*2*pi));
                % A_h the single harm matrix that relates laser to transd:
                A_h = Ld*inv(obj.Zred)*fT;
                ind = h+[0,1,2,3]*obj.Hmax;
                obj.Ahglobal(ind,ind) = A_h;
                if h==1
                    obj.Rtarget_h1 = Ld*inv(obj.Zred)*Ex_h1;
                end
            end
            obj = AssembleRtarget(obj,Rtarget_h1);
        end
        function R   = ResidualOuterLoopOld(obj,X)
            % Residual outer loop style
            % [LaserFS; LaserBS; LaserFC; LaserBC] = 
            %          Ahglobal*[TransFS; TransBS; TransFC; TransBC] + Rtarget;            
            ImposeVoltage(obj,X);           
            R=[...
                obj.rtc.par.in5_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_ave(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_ave(8+2:8+obj.Hmax+1)';...
                ] - ...
                obj.Ahglobal*[...
                obj.rtc.par.in3_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_ave(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in3_coeffs_ave(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_ave(8+2:8+obj.Hmax+1)';...
                ] + ...
                obj.Rtarget;
            noise_R=[...
                obj.rtc.par.in5_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in5_coeffs_var(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in6_coeffs_var(8+2:8+obj.Hmax+1)';...
                ] - ...
                obj.Ahglobal*[...
                obj.rtc.par.in3_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_var(0+2:0+obj.Hmax+1)';...
                obj.rtc.par.in3_coeffs_var(8+2:8+obj.Hmax+1)';...
                obj.rtc.par.in4_coeffs_var(8+2:8+obj.Hmax+1)';...
                ];
            disp('std(noise of R):')
            disp(sqrt(norm(noise_R)))
            if obj.scaling_residual
                R=R.*obj.scaling_residual;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end