%% Reads the value of the line corresponding to the name parameter
% After loading an ini file with load_ini_file(), you can extract any
% parameter using it parameter name. The value located after '=' will be
% returned.
% If the line does not exist, value is set to nan, unless a default value
% is provided.
% If an entry is found, but there is no value after '=', an error is raised
% If there is more than one entry, an error is raised unless you sepecify 
% and extra flag preceding this entry (in which case, the first entry after 
% that is returned)
%
% -------------------------------------------------------------------------
% Model : 
%   value = read_ini_value(ini_content, name, default_value, class_label)
%
% -------------------------------------------------------------------------
% Inputs : 
%   ini_content(Nx1 Cell Array) : a cell array returned by load_ini_file().
%
%   name(STR) : the name of the parameter to be read. T a partial name can
%               be used if there is no multiple entries with that name
%
%   default_value(any value) - optional - default is NaN : If no entry is 
%               found, default_value is returned
%
%   section_name(STR) - optional : if multiple entries are found, a search 
%               for section_name string is done, and the first valid entry
%               following section_name will be returned. To prevent errors,
%               use the full name of the class name, with the brackets.
%
%   separator(STR) - optional - default is '|' - If the type of data
%               imported is a list, split the results using this separator.
%               A list (even of numbers) must be surrounded by double
%               quotes.
%
% -------------------------------------------------------------------------
% Outputs :
%   value(BOOL, STR, FLOAT or N x 1 arrays) : If found, the value read 
%               after '=' is read otherwise default value is returned
%
% -------------------------------------------------------------------------
% Extra Notes:
%   If you want to import an array, your array must be surrounded by double
%   quotes. The default seprator is '|' but can be modified.
%   Final trailing whitespaces are ignored.
%   If a value ends with a semicolon, the character will be deleted during 
%   import in load_ini_file.m.
%
% -------------------------------------------------------------------------
% Author(s):
%   Antoine Valera  

function [value,exact_line] = read_ini_value(ini_content, name, default_value, section_name, separator)
   
    if nargin < 3 || isempty(default_value);default_value = NaN;        end
    if nargin < 4 || isempty(section_name); section_name = [];          end
    if nargin < 5 || isempty(separator);    separator = '|';            end

    %% Find lines matching 'name'
    Line = contains(ini_content,name);
    exact_line = [];
    
    %% Checking corner cases
    % -- multiple entries, but no rules for knowing which one to use
    if sum(Line) > 1 && (nargin < 4 || isempty(section_name))
        error_box(['more than one value found for ',name,' in ini file,',...
                   'but the [section] was not specified'], 0)
    end
    
    % -- line doesn't exist, and default_value is set
    if sum(Line) == 0 
        if nargin < 3
            fprintf(['no entry found for ', name,' in ini file ; setting',...
                     'value as NaN\n']);
        end
        value = default_value;
    end
    
    % -- multiple entries, and the appropriate class was indicated
    if sum(Line) > 1 && nargin == 4 && ~isempty(section_name)
       Header_location = find(contains(ini_content,section_name) == 1);
       if isempty(Header_location)
            error_box(['The specified section ',section_name, ' does not exist. Check for typo'], 0)
       end
       Line(1:Header_location) = 0; %zero any previous valid lines
       right_location = find(Line,1,'first'); %find index of right header
       Line(right_location+1:end) = 0; %zero any later valid lines
    end

    %% Now that all corner cases were solved, only one valid line is left
    if sum(Line) == 1 
        exact_line = find(Line == 1);
        txt = [ini_content{Line}];
        line_text = txt(strfind(txt, '= ')+1:end);
        if strcmpi(strtrim(line_text),'false') 
            value = false;
        elseif strcmpi(strtrim(line_text),'true')
            value = true;
        elseif isnumeric(str2double(line_text)) && ~isnan(str2double(line_text))
            value = str2double(line_text);
        elseif ~isempty(line_text)
            value = strip(strrep(strrep(line_text,'"',''),'''',''));
            if contains(value, separator)
                value = split(value,separator);
            end
        else
%             error_box(['entry value not identified for ', name,' in .ini file. There may be a missing parameter or a typo. Please check the ini file'], 0)
        end

    elseif sum(Line) == 0 && isnan(default_value)
%          error_box(['entry value not identified for ', name,' in .ini file. There may be a missing parameter or a typo, and no default value is provided'], 1)
    end
end