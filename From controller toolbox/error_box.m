%% Popup error/information message
% Popup window displaying a message. It can be used to catch an error or
% for any message. Chose a value for carry on to sue different symbols
%
% -------------------------------------------------------------------------
% Model : 
%   error_box(error_msg, carry_on)
%
% -------------------------------------------------------------------------
% Inputs : 
%   message(STR) : 
%                                   The message to display in the popup
%                                   window
%
%   carry_on(INT) - optional - default is true:  
%                                   If 0, raises an error message and
%                                   interrupt function
%                                   If 1, messge is a warning.
%                                   If 2, messge is a help dialog.
%                                   If 3, messge is a simple popup message.
%                                   otherwise the function will resume
%                                   after the popup closure.
%
% -------------------------------------------------------------------------
% Outputs :
% -------------------------------------------------------------------------
% Extra Notes:
% -------------------------------------------------------------------------
% Author(s):
%   Antoine Valera
% -------------------------------------------------------------------------
% Revision Date:
%   19-05-2018


function error_box(message, carry_on)
    if nargin < 2
        carry_on = 2;
    end
    
    %% Prompt message
    if ~carry_on
        uiwait(errordlg(message));   
    elseif carry_on == 1
        uiwait(warndlg(message)); 
    elseif carry_on == 2
        uiwait(helpdlg(message)); 
    elseif carry_on == 3
        uiwait(msgbox(message)); 
    end
    
    drawnow; pause(0.1);
    
    %% If carry_on is 0, interrupt running processes
    if ~carry_on
        error([message,'\n Process aborted']);
    end
end

