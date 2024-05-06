function fn_get_attenuation_exp()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all

%Reduce resolution of data
time_step = 1/50e6;
max_time = 3.5e-6;

%Load response variables
load("g4-s8_x-8000_z0_025.mat","exp_data");

figure;

total_transducer_els = 128;
for aperture_size = [4 8 16 32]
    aperture_vec = 1:aperture_size:total_transducer_els;

    %Loop over apertures
    for i = 1:(length(aperture_vec) - 1)

        %Define aperture
        aperture = aperture_vec(i):(aperture_vec(i + 1) - 1);
    
        %Calculate displacements
        aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
        aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);
        
        %Find response start
        ratio_of_max = 0.007;
        response_start_vec = find(abs(aperture_dsp_data) > ratio_of_max*max(abs(aperture_dsp_data)));
        response_start_vec_max = find(abs(aperture_dsp_data) > 1.2*ratio_of_max*max(abs(aperture_dsp_data)));
        for ii = 1:length(response_start_vec)
            if response_start_vec(ii) - response_start_vec_max(ii) < 1
                response_start = response_start_vec(ii);
                break
            end
        end

        translate_time = exp_data.time(response_start);
        exp_time_translated = exp_data.time - translate_time;
        %Calculate required indices
        half_time_step = 1/2 * time_step;
        exp_i = find(exp_time_translated > -half_time_step & exp_time_translated < (max_time + half_time_step));
    
        %Redefine time and displacements
        time_data = exp_time_translated(exp_i)';
        dsp_data = aperture_dsp_data(exp_i)';

        %Scale sim data
        dsp_norm = dsp_data/max(abs(dsp_data)); %Scale exp response to match sim response
        %Plot final results
        n_pts_half = round(length(dsp_norm)/2);
        attenuation = 1/max(dsp_norm(n_pts_half:end));
        attenuation_coef = 20*log10(attenuation) / 8;

        
        % porosity = (attenuation_coef - 0.5)/(5/6);
        
        % porosity = (attenuation_coef - 1.935854)/(0.139759);
        % fprintf("attenuation_coef: %.2f, porosity: %.2f\n", attenuation_coef, porosity)

        % plot(time_data, dsp_data, 'b')
        % hold on
    end
end





end

