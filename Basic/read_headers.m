function params = read_headers( filename )
% Read experimentHeader.ini - as generated using Labview for acquisition

params = [];
CurrSection = 'More';    
if isempty(filename)
    filename = uigetfile('', 'Choose experiment header file');
end

f = fopen(filename,'r');                % open file

while ~feof(f) 
    %Read next line (removing any white space
    newLine = strtrim(fgetl(f)); 
    if isempty( newLine ), continue;    end
    if newLine(1) == ';' || newLine(1) == '#',  continue;    end    % This is a comment line
    if ( newLine(1)=='[' ) && (newLine(end)==']' )                  % This is a new section
        CurrSection = matlab.lang.makeValidName( newLine(2:end-1) );
        params.(CurrSection) = [];
        if strcmpi(CurrSection, 'GLOBALPARAMETERS') || strcmpi(CurrSection, 'STATISTICS')
            makeNum = true;   
        else, makeNum = false; 
        end
    else                                         
        [fieldname,value] = strtok(newLine, '=');                   % This is a new field
        [fieldname, value] = check_params( fieldname, value, makeNum);
        if ~isempty(fieldname)
            params.(CurrSection).(fieldname) = value;
        end
    end
end

fclose(f);
end

function [ fld, val ] = check_params( fld, val, makeNum)
% check the integrity, rename, convert to numeric etc
fld = matlab.lang.makeValidName(fld);
eql = strfind(val, '=');
val = strtrim(val(eql+1:end));
if makeNum, val = str2num(val); end
end