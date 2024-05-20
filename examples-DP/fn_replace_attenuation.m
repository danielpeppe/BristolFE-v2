function fn_replace_attenuation()

% Specify the CSV file name
csvFileName = 'data/porosity_response_data_6.csv';

% Specify the column name or index in the CSV file containing the folder paths
targetColumn1 = 'batchID'; % You can change this to an index, e.g., 2
targetColumn2 = 'responseID';

% Read the CSV file into MATLAB
opts = detectImportOptions(csvFileName);
data = readtable(csvFileName, opts);

% Find the index of the target column if using column name
if ischar(targetColumn1)
    targetColumn1 = find(strcmp(opts.VariableNames, targetColumn1));
end

% Loop over each row in the target column
for i = 2:height(data)
    batch_num = data{i, targetColumn1}{1};
    response_num = data{i, targetColumn2}{1};
    
    mat_path = "data/" + string(batch_num) + "/" + string(response_num) + ".mat";
    
    % Check if the .MAT file exists
    if isfile(mat_path)
        % Load the .MAT file
        load(mat_path, 'op_config');
        o

        % Assume 'op_config' contains a variable named 'value' we need
        if isfield(op_config, 'value')
            % Update the CSV table with the value from .MAT file
            data.value(i) = op_config.value;
        else
            disp(['"value" not found in ', matFilePath]);
        end
    else
        disp([matFilePath, ' does not exist']);
    end
end

% Write the updated table back to the CSV file
writetable(data, 'updated_output.csv');
disp('Updated CSV file has been saved as "updated_output.csv"');


end