classdef it_2out_num_interface < rtc_interface
    %NTMD_INTERFACE  Interface to the real-time control device that controls
    %  the nonlinear tuned mass damper experiment.
    
    % Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015
    
    properties
        fourier;
        datafields;
    end
    
    methods
        function obj = it_2out_num_interface()
            %NTMD_INTERFACE  Construct an NTMD_INTERFACE object.
            
            % Indices into the array of Fourier variables
            n_coeff = length(obj.par.in5_coeffs);
            obj.fourier.n_modes = n_coeff/2 - 1;
            obj.fourier.n_ave = obj.par_info(obj.par_idx('in5_coeffs_arr')).count/n_coeff; % Number of periods that are averaged to get the result
            obj.fourier.idx_DC = n_coeff/2 + 1;
            obj.fourier.idx_AC = setdiff(2:n_coeff, obj.fourier.idx_DC);
            obj.fourier.idx_fund = [2, obj.fourier.idx_DC + 1];
            obj.fourier.idx_higher = [obj.fourier.idx_fund(1) + (1:obj.fourier.n_modes - 1), ...
                                      obj.fourier.idx_fund(2) + (1:obj.fourier.n_modes - 1)];
            obj.fourier.idx_sin = 1 + (1:obj.fourier.n_modes);
            obj.fourier.idx_cos = n_coeff/2 + 1 + (1:obj.fourier.n_modes);
            
            % Default options for the experiment
            obj.opt.samples = 100; % Number of samples to record
            obj.opt.downsample = 0; % Number of samples to ignore for every sample recorded
            obj.opt.wait_time = 1; % Time (in secs) to wait for Fourier coefficients to settle
            obj.opt.max_waits = 10; % Maximum number of times to wait
            obj.opt.max_picard_iter = 5; % Maximum number of Picard iterations to do
            obj.opt.x_coeffs_var_tol_rel = 1e-3; % Maximum (normalised) variance of Fourier coefficients for steady-state behaviour
            obj.opt.x_coeffs_var_tol_abs = 2e-3; % Maximum (absolute) variance of Fourier coefficients for steady-state behaviour
            obj.opt.x_coeffs_tol = 1e-1; % Maximum tolerance for difference between two Fourier coefficients (mm)
            
            % Data recording fields
            obj.datafields.stream_id = 1; % The stream to use for data recording
            obj.datafields.static_fields = {'input_filter_freq',...
                                'sample_freq', 'firmware'};
            obj.datafields.dynamic_fields = {'voltage_freq',...
                                'voltage7_coeffs','voltage8_coeffs', ...
                                'in5_coeffs_ave', 'in5_coeffs_var',...
                                'in6_coeffs_ave', 'in6_coeffs_var', ...
                                'in1_coeffs_ave', 'in1_coeffs_var',...
                                'in2_coeffs_ave', 'in2_coeffs_var', ...                             
                                'in3_coeffs_ave', 'in3_coeffs_var',...
                                'in4_coeffs_ave', 'in4_coeffs_var'};
            obj.datafields.stream_fields = {'time_mod_2pi', ...
                                'in5', 'in6', ...
                                'in1', 'in2', ...
                                'in3', 'in4','voltage7','voltage8'}; 
            
            % Set the filter to 50 Hz
            obj.par.input_filter_freq = 50/obj.par.sample_freq;
            obj.par.rand_filter_freq = 10/obj.par.sample_freq;
        end
    end
    
end

