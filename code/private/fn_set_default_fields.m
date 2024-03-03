function new_struct = fn_set_default_fields(old_struct, default_struct, varargin)
%USAGE
%	new_struct = fn_set_default_fields(old_struct, default_struct);
%SUMMARY
%	Use to add default fields and values to a structured variable, such as
%	options for a function.
%AUTHOR
%	Paul Wilcox (Dec 2003)
%INPUTS
%	old_struct - original structured variable
%	default_struct - structured variable containing default fields and
%	values
%OUTPUTS
%	new_struct - updated structured variable. All existing fields and values
%	in old_struct will be preserved, but any fields found in default_struct
%	and their values will be added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
	using_for_op = [];
else
	using_for_op = varargin{1};
end

new_struct = old_struct;
default_fieldnames = fieldnames(default_struct);
for ii=1:length(default_fieldnames)
	if ~isfield(new_struct, default_fieldnames{ii})
		new_struct = setfield(new_struct, default_fieldnames{ii}, getfield(default_struct, default_fieldnames{ii}));
    end
end

%Checks for options
new_struct_fieldnames = fieldnames(new_struct);
if using_for_op
    for ii=1:length(new_struct_fieldnames)
        %Check if default is defined
        if ~isfield(default_struct, new_struct_fieldnames{ii})
            fprintf("CAUTION: default_op not defined for %s\n", new_struct_fieldnames{ii})
        else
            %Print non-default field values
            default_value = default_struct.(new_struct_fieldnames{ii});
            new_struct_value = new_struct.(new_struct_fieldnames{ii});
            % for numeric values
            if isnumeric(new_struct_value)
                if default_value ~= new_struct_value
                    fprintf('op.%s: %.4f (default: %.4f)\n', new_struct_fieldnames{ii}, new_struct_value, default_value);
                end
            
            elseif ischar(new_struct_value) || isstring(new_struct_value)
                if ~strcmpi(default_value, new_struct_value)
                    fprintf('op.%s: %s (default: %s)\n', new_struct_fieldnames{ii}, new_struct_value, default_value);
                end
            % % for matrices TBC
            % elseif ismatrix(new_struct_value) && ~(ischar(new_struct_value) || isstring(new_struct_value))
            %     if default_value ~= new_struct_value
            %         fprintf('op.%s: (non-default)\n', new_struct_fieldnames{ii}); %hard to deal with matrices
            %     end
            % 
            else
                error('Unsupported variable type in option');
            end
        end
    end
    %%Enable the following when required
    % %Finally check for any options set in default that aren't present in
    % %old_struct
    % for ii=1:length(default_fieldnames)
    %     %Check if new_struct is defined
    %     if ~isfield(old_struct, default_fieldnames{ii})
    %         fprintf("CAUTION: %s options is not defined in default options\n", default_fieldnames{ii})
    %     end
    % end
end
return
