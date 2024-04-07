function fn_process_mat_data(batch_num_lower, varargin)

%Define upper batch bound as lower batch bound if nothing input
if isempty(varargin)
    batch_num_upper = batch_num_lower;
else
    batch_num_upper = varargin{1};
end

%Create csv file for porosity_response_data
data_path = "data/";
file_ID = fopen(data_path + "porosity_response_data.csv", 'w');

%Reduce resolution of data
time_step = 1/50e6;
max_time = 3.5e-6;
new_time_data = 0:time_step:max_time;

%Write time data to top of csv
fprintf(file_ID, '%s,%s,%s,', 'batchID', 'responseID', 'porosity');
for ii = 1:(length(new_time_data) - 1)
    fprintf(file_ID, '%s,', "timeData" + ii);
end
fprintf(file_ID, '\n');
fprintf(file_ID, '%s,%s,%s,', 'NaN', 'NaN', 'NaN');
for ii = 1:(length(new_time_data) - 1)
    fprintf(file_ID, '%.15f,', new_time_data(ii));
end
fprintf(file_ID, '\n');

%Loop over batches
for batch_number = batch_num_lower:batch_num_upper
    %Get batch folder path
    folder_path = data_path + "batch" + batch_number;

    %Count the number of responses in batch
    responses = dir(fullfile(folder_path, '*.mat'));
    n_responses = numel(responses);
    fprintf('writing batch: %d with responses: %d...', batch_number, n_responses)

    %Loop over responses
    for i = 1:n_responses
        %Load response variables
        file_path = folder_path + "/response" + i + ".mat";
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
        fprintf(file_ID, '%f,', op_config.porosity);
        for ii = 1:length(new_dsp_data)
            fprintf(file_ID, '%.15f,', new_dsp_data(ii));
        end
        fprintf(file_ID, '\n');
    end
    fprintf('completed\n')
end

fclose(file_ID);

end

