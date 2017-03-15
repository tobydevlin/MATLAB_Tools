% /////// checkme ///////
%
% function checkme(S,resfil,variables_arguments)
%
% function to determine whether fv_get_dave or its sister fv_get_... functions are being used correctly.
% The danger in using fv_get_dave and its sister functions is that you may call the output
% stucture a 2nd, 3rd .. nth time with different inputs but with the same
% variable name for the output structure. If this happens fv_get_dave
% ignores the inputs and works on the information it has saved in its
% previous calling!

function checkme(S,fname,resfil,variable_arguments)

if ischar(resfil)
    if ~strcmp(resfil,S.resfil)
        error([fname ' has stored information from a previous call which you do not want to use - you must clear your output structure'])
    end
else
    if ~isequal(resfil,S.resfil)
        error([fname ' has stored information from a previous call which you do not want to use - you must clear your output structure'])
    end
end

if ~isequal(variable_arguments,S.varargin)
    % -- data
    [i1, l1] = ismember('data',variable_arguments(1:2:end));
    [i2, l2] = ismember('data',S.varargin(1:2:end));
    if all([i1 i2])
        if size(variable_arguments{2*l1}) ~= size(S.varargin{2*l2})
             error([fname ' has stored information from a previous call which you do not want to use - you must clear your output structure'])
        end
    elseif any([i1 i2])
        error([fname ' has stored information from a previous call which you do not want to use - you must clear your output structure'])
    end
    
    % -- stat
    [i1, l1] = ismember('stat',variable_arguments(1:2:end));
    [i2, l2] = ismember('stat',S.varargin(1:2:end));
    if all([i1 i2])
        if size(variable_arguments{2*l1}) ~= size(S.varargin{2*l2})
           error([fname ' has stored information from a previous call which you do not want to use - you must clear your output structure']) 
        end
   elseif any([i1 i2])
        error([fname ' has stored information from a previous call which you do not want to use - you must clear your output structure'])
    end
end
