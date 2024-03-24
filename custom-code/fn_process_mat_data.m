function fn_process_mat_data(batch_name)

% Prepare the data row

% Specify the file name
file_name = 'data.csv';

% Open the file in write mode
file_id = fopen(file_name, 'w');

% Write the data row to the file
% The format '%f;' specifies that each number should be written as a floating point,
% followed by a semicolon. The '\n' specifies a new line at the end of the row.
fprintf(file_id, '%f;', data_row);

% Optionally, you might want to write a new line at the end if you plan to add more rows later
fprintf(file_id, '\n');

% Close the file
fclose(file_id);

% Note: This will write the data row to 'data.csv' with semicolon delimiters.
% If you want to include column headers, you can write them to the file before writing the data.


folder_path = ['data/' batch_name];

% Count the number of .mat files
mat_files = dir(fullfile(folder_path, '*.mat'));
n_mat_files = numel(mat_files);

% New sampling rate is 1/n of the original
n = 1;
[p, q] = rat(1/n); % For a downsampling factor n

for i = 1:n_mat_files
    % Step 1: Load the .mat file
    % Replace 'yourFile.mat' with the path to your .mat file
    load([folder_path '/response' char(string(i)) '.mat'], 'res', 'steps', 'op_save')
    
    %Get time and displacement data
    time_data = steps{1,1}{1}.load.time;
    dsp_sum = sum(res{1,1}{1}.dsps);
    dsp_norm = dsp_sum/max(abs(dsp_sum));

    %Downsample results
    time_data_ds = resample(time_data, p, q);
    dsp_data_ds = resample(dsp_norm, p, q);

    %Get porosity level
    porosity = op.save{1,1}{1}.porosity;

    data_row = [porosity, time_data]; % Concatenate porosity and time data

    writematrix(timeData, ['data_' bach_name '.csv']);

end
end

