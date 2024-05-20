function fn_process_mat_data_exp()

close all

%Create csv file for porosity_response_data
file_ID = fopen("exp_response_data.csv", 'w');

%Reduce resolution of data
time_step = 1/50e6;
max_time = 3.5e-6;
new_time_data = 0:time_step:max_time;

%Write time data to top of csv
fprintf(file_ID, '%s,%s,%s,%s,', 'batchID', 'responseID', 'porosity', 'attenuation-coef_dB-mm');
for ii = 1:length(new_time_data)
    fprintf(file_ID, '%s,', "timeData" + ii);
end
fprintf(file_ID, '\n');
fprintf(file_ID, '%s,%s,%s,%s,', 'NaN', 'NaN', 'NaN', 'NaN');
for ii = 1:length(new_time_data)
    fprintf(file_ID, '%.15f,', new_time_data(ii));
end
fprintf(file_ID, '\n');

%Load response variables
load("g4-s8_x-8000_z0_025.mat","exp_data");

figure;

total_transducer_els = 128;
for aperture_size = [4 8 16 32]
    aperture_vec = 1:aperture_size:total_transducer_els;

    %Loop over apertures
    for i = 1:(length(aperture_vec) - 1)
    
        %Write response ID
        fprintf(file_ID, '%s,%s,', "exp_aperture" + aperture_size, "response" + i);
        %Porosity
        fprintf(file_ID, '%f,', 0);
    
        %Define aperture
        aperture = aperture_vec(i):(aperture_vec(i + 1) - 1);
    
        %Calculate displacements
        aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
        aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);
        
        %Find response start
        ratio_of_max = 0.007;
        % if aperture_size == 8 || aperture_size == 16
        %     % ratio_of_max = 0.0079;
        % elseif aperture_size == 32
        %     % ratio_of_max = 0.0020;
        %     disp('hello')
        % end
        response_start_vec = find(abs(aperture_dsp_data) > ratio_of_max*max(abs(aperture_dsp_data)));
        response_start_vec_max = find(abs(aperture_dsp_data) > 1.2*ratio_of_max*max(abs(aperture_dsp_data)));
        for ii = 1:length(response_start_vec)
            if response_start_vec(ii) - response_start_vec_max(ii) < 1
                response_start = response_start_vec(ii);
                break
            end
        end
        %response_start_vec_cleaned = response_start_vec(diff(response_start_vec) < 2);
        %response_start = response_start_vec_cleaned(1);
        %respose_start_vec = response_start_vec(response_start_vec > 50);
        %response_start = respose_start_vec(1);
        translate_time = exp_data.time(response_start);
        exp_time_translated = exp_data.time - translate_time;
        %Calculate required indices
        half_time_step = 1/2 * time_step;
        exp_i = find(exp_time_translated > -half_time_step & exp_time_translated < (max_time + half_time_step));
    
        %Redefine time and displacements
        time_data = exp_time_translated(exp_i)';
        dsp_data = aperture_dsp_data(exp_i)';

        % %Reduce resolution of displacement data
        % new_dsp_data = interp1(time_data, dsp_data, new_time_data, 'linear');
        % %Cleaning
        % new_dsp_data(1) = 0;
        % new_dsp_data = new_dsp_data(1:(end-1)); %last value is NaN
        % new_time_data = new_time_data(1:(end-1));

        %Scale sim data
        dsp_norm = dsp_data/max(abs(dsp_data)); %Scale exp response to match sim response
        %Plot final results
        n_pts_half = round(length(dsp_norm)/2);
        attenuation = 1/max(dsp_norm(n_pts_half:end));
        attenuation_coef = 20*log10(attenuation) / 8;
        fprintf(file_ID, '%s,', attenuation_coef);


        %Write Displacements
        for ii = 1:length(dsp_data)
            fprintf(file_ID, '%.15f,', dsp_data(ii));
        end
        fprintf(file_ID, '\n');
            
        plot(time_data, dsp_data, 'b')
        hold on
        % plot(new_time_data, dsp_data, 'r')
        % hold on
    end
end
fclose(file_ID);

end
