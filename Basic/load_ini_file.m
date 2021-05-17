%% Read ini file OR file paths if 'section' is indicated 
% You can either directly read the content of an ini file by passing the
% full adress of the file, or you can pass an ini file that contains a list
% of paths, plus the section name that precede this path.
% see details about ini stucture in the Extra Notes section
%
%
% -------------------------------------------------------------------------
% Model : 
%   ini_info = load_ini_file(ini_file_path, section)
%
% -------------------------------------------------------------------------
% Inputs : 
%   ini_file_path(STR) - optional : the full path of an ini file ;
%               OR the path of an ini file that contains pass, in which
%               case you need to pass the section parameter too
%               OR [], in which case if an existing configuration_file_path  
%               .ini is present in the ./Core folder, it will be used
%
%   section(STR) - optional : if you pass a file with paths in it, you must
%               indicate the name of the section that contains your path of
%               interest
%
% -------------------------------------------------------------------------
% Outputs :
%   ini_info(Nx1 Cell array) : Cell array axtracted from the ini file.
%               There is one cell per non-empty line
%
% -------------------------------------------------------------------------
% Extra Notes:
%   #### Case 1
% 	full .ini file structure must be in the form
%
%     [Section name with numbers]
%     Variable 1 name = variable 1 numeric value  
%     Variable 2 name = variable 2 boolean value
% 
%     [Section name with strings]
%     Variable 3 name = "variable 3 whatever string" 
%     Variable 4 name = 'variable 4 whatever string' 
%     Variable 5 name = variable 5 whatever string 
%
%   #### Case 2
% 	path .ini file structure must be in the form
%
%     [Section name for path 1]
%     X:\Wherever that file is\it_is_there.ini
%     [Section name for path 2]
%     Y:\Wherever that file is\it_is_another_one.ini
%
% -------------------------------------------------------------------------
% Author(s):
%   Antoine Valera

function [ini_content, filepath] = load_ini_file(ini_file_path, section)
    %% identify filename
    if nargin < 1 || isempty(ini_file_path)
        ini_file_path = '@/configuration_file_path.ini';
    end
    if nargin < 2 || isempty(section)
        filepath = ini_file_path;
    else
        [filepath, ini_file_path] = get_setup_ini_path(ini_file_path, section);
        if isempty(filepath)
            error_box([ini_file_path(2:end), ' file do not exist and was created. Please restart the controller'], 0)
        elseif filepath == -1
            error_box(sprintf('The Section "%s" in the ini file "%s" could not be found, Check that the .ini file is correct and points to the right destination.',section,ini_file_path), 0)
        end
    end
    
    %% open ini file ; read values
    fid = fopen(filepath,'r');
    if fid == -1
        if ~exist(filepath) && exist(strrep(filepath,'.ini','.template'))
            copyfile(strrep(filepath,'.ini','.template'),filepath);
            error_box([ini_file_path(2:end), ' file do not exist and was created from the existing .template file. Please check the ini file content to adjust it to the current setup and restart the controller'], 1);
            fid = fopen(filepath,'r');
        else
            error_box('Configuration file not found. The path indicated in the ini file is probably wrong', 0);
        end
    end
    text = textscan(fid, '%s','Delimiter','');
    fclose(fid);
    text = text{1};
    
    %% split lines
    %ini_content = text{:};
    ini_content = cleanup_text(text); 
end

function [config_file_path, ini_file_path] = get_setup_ini_path(ini_file_path, section)
    %% Get the line following section if you passed an ini files with paths
    
    if strcmp(ini_file_path(1),'@')    
        % This requires an configuration_file_path.ini file specifically in
        % the ./Core folder. It's not a good generic solution
        [~, ~, ~, ~, ini_folder_path] = adjust_pathnames(mfilename('fullpath'));
        ini_file_path = [ini_folder_path, ini_file_path(3:end)];
    else
        %% Correct paths for typo
        [ini_file_path, ~, ~, ~, ~] = adjust_pathnames(ini_file_path);
    end
    
    %% Open path file
    path_file = fopen(ini_file_path);
    if path_file == -1
        %% If it doesn't exist, it must be generated
        create_initial_file(ini_file_path);
        path_file = fopen(ini_file_path);
    end
    
    %% Find the right line, with the path
    config_file_path = find_right_header(path_file, section); 
    fclose(path_file);
end

function result = find_right_header(path_file, match)
    %% identify a line following indicated header
    tline = 'whatever';
    while ischar(tline) && ~strcmp(tline, match)
        tline = fgetl(path_file);
    end
    result = fgetl(path_file);
end

function text = cleanup_text(text)
    %% normalize text to remove trailing zero and double whitespaces
    for idx = 1:size(text,1)
        line = text{idx};

        % removes heading and leading zeroes
        line = strip(line);

        % remove double spaces
        while contains(line,'  ')
            line = strrep(line,'  ',' ');
        end

        % make sure there is one space before and after =
        if ~contains(line,'= ')
            line = strrep(line,'=','= ');
        end
        if ~contains(line,' =')
            line = strrep(line,'=',' =');
        end

        semicolon_idx = strfind(line,';');
        if ~isempty(semicolon_idx) && semicolon_idx == size(line,2)
            line = line(1:end-1);
        end

        text{idx} = line;
    end
end

function create_initial_file(ini_file_path)

    error_box( [ini_file_path,...
               ' was not found. If it is your first use on that computer ',...
               ' or if you just reinstalled the controller, you have to update the content of this .ini file.',...
               ' The file will now be created from a standard .template file but you have to check the .ini content',...
               ' Please open the file and adjust the relevant variables and paths.',...
               ' - For analysis mode (offline), you just have to update the vaa3d path.',...
               ' - For acquisition mode (online), Please calibrate the relevant objectives, adjust communications ports',...
               ' for all required devices, the ethernet adresses. Check AOL configuration, in particular the max power passed to crystals'], 1)   
    if ~exist(ini_file_path, 'file') && exist(strrep(ini_file_path,'.ini','.template'), 'file')
        copyfile(strrep(ini_file_path,'.ini','.template'),ini_file_path);
        fid = fopen(ini_file_path,'r');
        t = textscan(fid, '%s','Delimiter',''); t = t{1};
        p = mfilename('fullpath');
        idx = sort([strfind(p,'/'),strfind(p,'\')]);
        p1 = [p(1:idx(end-1)),'+default\setup.ini'];
        p2 = [p(1:idx(end-1)),'+default\calibration.ini'];
        t{2} = p1;
        t{4} = p2;
        fclose(fid);
        fid = fopen(ini_file_path,'wt');
        fprintf(fid,'%s\n',t{:});
        fclose(fid);
    end
end

