function fn_print_default_options(op, default_op)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%turn off warning backtrace
warning('OFF', 'BACKTRACE');

new_struct = op;
default_struct = default_op;
new_struct_fieldnames = fieldnames(new_struct);
for ii=1:length(new_struct_fieldnames)
    %Check if default is defined
    if ~isfield(default_struct, new_struct_fieldnames{ii})
        warning("default_op not defined for %s", new_struct_fieldnames{ii})
    else
        %Print non-default field values
        default_value = default_struct.(new_struct_fieldnames{ii});
        new_struct_value = new_struct.(new_struct_fieldnames{ii});
        % for parameters for multiple sims
        if strcmpi(new_struct_fieldnames{ii},'params')
            %do nothing
        % for numeric values
        elseif isnumeric(new_struct_value)
            if default_value ~= new_struct_value
                fprintf('op.%s: %.4f (default: %.4f)\n', new_struct_fieldnames{ii}, new_struct_value, default_value);
            end
        elseif ischar(new_struct_value) || isstring(new_struct_value)
            if ~strcmpi(default_value, new_struct_value)
                fprintf('op.%s: %s (default: %s)\n', new_struct_fieldnames{ii}, new_struct_value, default_value);
            end
        elseif iscell(new_struct_value)
            if ~isempty(new_struct_value)
                fprintf('op.%s: custom values set\n', new_struct_fieldnames{ii});
            end
        else
            error('Unsupported variable type in option');
        end
    end
end


end

