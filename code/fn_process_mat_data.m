function fn_process_mat_data(varargin)

just_count_responses = 1;

%Define upper batch bound as lower batch bound if nothing input
if isempty(varargin)
    batch_num_lower = 1;
    batch_num_upper = 1000;
else
    batch_num_lower = varargin{1};
    batch_num_upper = varargin{2};
end

%Create csv file for porosity_response_data
data_path = "../data/";
file_ID = fopen(data_path + "porosity_response_data.csv", 'w');

%Reduce resolution of data
time_step = 1/50e6;
max_time = 3.5e-6;
new_time_data = 0:time_step:max_time;

%Write time data to top of csv
fprintf(file_ID, '%s,%s,%s,%s,', 'batchID', 'responseID', 'porosity', 'attenuation-coef_dB-mm');
for ii = 1:(length(new_time_data) - 1)
    fprintf(file_ID, '%s,', "timeData" + ii);
end
fprintf(file_ID, '\n');
fprintf(file_ID, '%s,%s,%s,%s,', 'NaN', 'NaN', 'NaN','NaN');
for ii = 1:(length(new_time_data) - 1)
    fprintf(file_ID, '%.15f,', new_time_data(ii));
end
fprintf(file_ID, '\n');

%Loop over batches
total_n_responses = 0;
for batch_number = batch_num_lower:batch_num_upper
    %Get batch folder path
    batch_path = data_path + "batch" + batch_number;
    if isfolder(batch_path)
        %Count the number of responses in batch
        responses = dir(fullfile(batch_path, '*.mat'));
        n_responses = numel(responses);
        total_n_responses = total_n_responses + n_responses;    
        %Loop over responses
        if ~just_count_responses
            fprintf('writing batch: %d with responses: %d...\n', batch_number, n_responses)
            for i = 1:n_responses
                %Load response variables
                file_path = batch_path + "/response" + i + ".mat";
                load(file_path, 'dsp_data', 'time_data', 'op_config')
                if size(dsp_data) ~= size(time_data)
                    warning('Time and Displacement data are different sizes')
                end
                
                %Reduce resolution of displacement data
                new_dsp_data = interp1(time_data, dsp_data, new_time_data, 'linear');
                new_dsp_data = new_dsp_data(1:(end-1)); %last value is NaN
        
                %Write batch ID
                fprintf(file_ID, '%s,%s,', "batch" + batch_number, "response" + i);
        
                %Write time data and porosity
                % fprintf(file_ID, '%s,', "time" + i);
                % fprintf(file_ID, '%f,', op_config.porosity);
                % for ii = 1:length(new_time_data)
                %     fprintf(file_ID, '%.15f,', new_time_data(ii));
                % end
                % fprintf(file_ID, '\n');
        
                %Write displacement data (ensure NaN is placed to line up porosity)
                % fprintf(file_ID, '%s,', 'NaN');
                % fprintf(file_ID, '%s,', 'NaN');
                % fprintf(file_ID, '%s,', "displacement" + i);
                % fprintf(file_ID, '%s,', 'NaN');
                
                %Write porosity
                fprintf(file_ID, '%f,', op_config.porosity);
                
                %Write attenuation
                fprintf(file_ID, '%f,', 20*log10(op_config.attenuation) / (2 * op_config.specimen_size * 1000));
    
                for ii = 1:length(new_dsp_data)
                    fprintf(file_ID, '%.15f,', new_dsp_data(ii));
                end
                fprintf(file_ID, '\n');
            end
        end
    end
end

fprintf('Total Number of responses: %d\n', total_n_responses)
fclose(file_ID);

end

