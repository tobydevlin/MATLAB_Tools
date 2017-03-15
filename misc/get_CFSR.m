% Get CFSR
% Submits Requests for CFSR and CFSv2 Data from RDA UCAR
% Requires a valid login which can be requested at:
%
% https://rda.ucar.edu/index.html?hash=data_user&action=register
%
%----------------------------------------
% Time periods available are 1979 - 2010 inclusive (CFSR)
% and 2011 - now (ish) (CFSv2)
%
% You can request both at the same time and it will submit two
% requests. But as they are on different grids you will need 
% to process the data separately (so put them in
% different folders and use the processing scripts).
%
%--------------------------
% VARIABLES--
%
% rain: Total Accumulation within an hour
% temp: Air temp
% rads: long & short wave averages 
%       over differing periods (0-1, 0-2, etc.)
%       needs to be adjusted to hourly averages
% rhum: relative humidity 
% wind: x and y component of wind
%
%--------------------------
% 
% INPUTS: (example)
%
% variables = {'wind';'mslp';'temp';'rhum';'rads';'rain'};
% ts = '01/12/2015';          % dd/mm/yyyy
% te = '01/12/2016';          % dd/mm/yyyy
% lon = [114 116];         		% min max
% lat = [-34 -30];         		% min max
%
% get_CFSR(variables, ts, te, lon, lat)
%-----------------------------------------
% -- TD Mar2015
%----------------------------------------
%
% Updates:
%
% 13/03/2017: updated to handle https on new NCAR RDA server (TDEVLIN)
%
%

function get_CFSR(variables, ts, te, lon, lat)


% -- DONT EDIT (Unless you know what to do)
%----------------------------------------
    
    % -- RDA codes
        code(1).rads.var = '3!7-0.2-1:0.5.192,3!7-0.2-1:0.4.192';
        code(1).rads.tim = '944,950,650,951,952,652';
        code(1).rads.grd = '62';
        code(1).rads.lvl = '107';

        code(1).rain.var = '3!7-0.2-1:0.1.8';
        code(1).rain.tim = '604,824,828,830,834,838';
        code(1).rain.grd = '57';
        code(1).rain.lvl = '107';

        code(1).wind.var = '3!7-0.2-1:0.2.2,3!7-0.2-1:0.2.3';
        code(1).wind.tim = '486,488,119,490,492,3';
        code(1).wind.grd = '62';
        code(1).wind.lvl = '223';

        code(1).temp.var = '3!7-0.2-1:0.0.0';
        code(1).temp.tim = '486,488,119,490,492,3';
        code(1).temp.grd = '62';
        code(1).temp.lvl = '107,210';

        code(1).mslp.var = '3!7-0.2-1:0.3.1';
        code(1).mslp.tim = '1,486,488,119,490,492';
        code(1).mslp.grd = '57';
        code(1).mslp.lvl = '219';

        code(1).rhum.var = '3!7-0.2-1:0.1.1';
        code(1).rhum.tim = '486,488,119,490,492,3';
        code(1).rhum.grd = '57';
        code(1).rhum.lvl = '221';


        code(2) = code(1);
        code(2).rads.grd = '68';
        code(2).temp.grd = '68';
        code(2).wind.grd = '68';

        ts = datenum(ts,'dd/mm/yyyy');
        te = datenum(te,'dd/mm/yyyy');

        lon = [floor(lon(1)) , ceil(lon(2))];
        lat = [floor(lat(1)) , ceil(lat(2))];

    if ~isunix
		    error('This requires the ''curl'' function (with openssl support), so needs to run on a linux machine')
		end	

    % -- Check Times
        
        if ts<datenum(2011,1,1)
            
            if te>datenum(2011,1,1)
                models = [1;2];
                warndlg('This is requesting data from two different models. They exist on different grids so must be processed separately')
            elseif ts<=datenum(2011,1,1)
                models = 1;
            end
            ts = max(ts,datenum(1979,1,1));
            
        elseif ts>=datenum(2011,1,1)
            
            models = 2;

        end

        dsid = {'ds093.1';'ds094.1'};
        ds   = {'093.1';'094.1'};

        % -- LOGIN and Request DATA
        answer{1} = '';
        answer{2} = '';
        [answer{2},answer{1}]=passwordEntryDialog('enterUserName',true,...
                                            'CheckPasswordLength',false);
try
    if exist('rda_auth_cookies','file')
        dat = importdata('rda_auth_cookies');
        if ~strcmpi(dat{end}(46:end),answer{1})
            
            system(['wget --save-cookies rda_auth_cookies --post-data="email='...
                    answer{1} '&passwd=' answer{2} ...
                    '&action=login" https://rda.ucar.edu/cgi-bin/login']);
                
        end
    else
        
        system(['wget --save-cookies rda_auth_cookies --post-data="email=' ...
                answer{1} '&passwd=' answer{2} ...
                '&action=login" https://rda.ucar.edu/cgi-bin/login']);
            
    end
catch 
    clearvars answer
    error('Couldnt get Cookies')
end


clearvars answer

for aa = models
    for bb = 1 : length(variables)

        if isunix
            name = strrep(code(aa).(variables{bb}).var,'!','%21');
        else
            name = code(aa).(variables{bb}).var;
        end
        time = code(aa).(variables{bb}).tim;
        grid = code(aa).(variables{bb}).grd;
        level = code(aa).(variables{bb}).lvl;

        % Hack to get matlab to ignore its own library path (matlab has an old libcurl w/o openssl on matlab 2013)
        % Set library path empty
        oldlib = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH','')
        
        % Setup https request
        sub = ['curl -b rda_auth_cookies -d "dsid=' dsid{aa}...
            '&rtype=S&rinfo=dsnum=' ds{aa} ...
            ';startdate=' datestr(ts,'yyyy-mm-dd HH:MM') ...
            ';enddate=' datestr(te,'yyyy-mm-dd HH:MM') ...
            ';product=' time ...
            ';parameters=' name ...
            ';level=' level ...
            ';grid_definition=' grid ...
            ';slat=' num2str(lat(1)) ...
            ';nlat=' num2str(lat(2)) ...
            ';wlon=' num2str(lon(1)) ...
            ';elon=' num2str(lon(2)) ...
            ';ofmt=netCDF" https://rda.ucar.edu/php/dsrqst.php'];

        system(sub,'-echo')
        
        % Set library path back (just in case its needed elsewhere)
        setenv('LD_LIBRARY_PATH',oldlib)
    end
    % Cleanup the cookies
    !rm login
    !rm rda_auth_cookies
end



function [Password, UserName] = passwordEntryDialog(varargin)
% PASSWORDENTRYDIALOG
% [Password, UserName] = passwordEntryDialog(varargin)
%
% Create a password entry dialog for entering a password that is visibly
% hidden.  Java must be enabled for this function to work properly.
%
% It has only been tested on the Windows platform in R2008a.  It should
% work in R2007a or later.
%
% The password box is created using the Java Swing component
% JPasswordField.
%
% Optional Input Arguments
% ------------------------
%
% 'enterUserName'       DEFAULT: false
% Display the user name entry box.  The user name entered must be at least
% one character or an error dialog is displayed.
%
% 'DefaultUserName'     DEFAULT: ''
% String value of user name to populate in User Name box upon creation.
%
% 'ValidatePassword'    DEFAULT: false
% Display dialog box to reenter password for validation purposes.
%
% 'CheckPasswordLength' DEFAULT: true
% Check the password length to ensure it meets the specified criteria.
%
% 'PasswordLengthMin'   DEFAULT: 2
% Minimum password length allowed.
%
% 'PasswordLengthMax'   DEFAULT: 8
% Maximum password length allowed.
%
% 'WindowName'          DEFAULT: 'Login'
% Title of the password entry window.
%
% Examples
% --------
%
% Create a password dialog box with the default options.
% -----------------------------------------------------------------------
% [Password] = passwordEntryDialog;
%
% Create a user name and password entry dialog box without password
% verification.
% -----------------------------------------------------------------------
% [Password, UserName] = passwordEntryDialog('enterUserName', true)
%
% Create a user name and password entry dialog box without password
% verification.  Set the user name to 'jdoe' upon startup.
% -----------------------------------------------------------------------
% [Password, UserName] = passwordEntryDialog('enterUserName', true,...
%      'DefaultUserName', 'jdoe')
%
% Create a password dialog box with password validation
% -----------------------------------------------------------------------
% [Password] = passwordEntryDialog('ValidatePassword', true);
%
% Create a user name and password entry dialog box with password
% verification.
% -----------------------------------------------------------------------
% [Password, UserName] = passwordEntryDialog('enterUserName', true,...
%      'ValidatePassword', true)
%
% Check the length of the password to be between 5 and 8 characters
% -----------------------------------------------------------------------
% [Password, UserName] = passwordEntryDialog('CheckPasswordLength', true,...
%      'PasswordLengthMin', 5,...
%      'PasswordLengthMax', 8)
%
% -----------------------------------------------------------------------
% Copyright (C) 2007-2008, Jesse B. Lai
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>

%% History
% The history of this program is outlined in this section.
%
% 20080616 - JBL - 0.0.2
% Started 20080612
%
% Updated to remove requirement for uicomponent function written by Yair
% Altman.  Now, the Java components are created manually and the
% undocumented feature javacomponent is used.
%
% A focus issue was worked out by using the drawnow command in a couple of
% places to allow the objects to be focused properly upon startup.
%
% 20080427 - JBL - 0.0.1
% Started 20080425.
%
% Initial version.  Uses Java to create the password entry box.  An edit
% box was initially attemped with the Java frame, but was occasionally
% getting Java exceptions when trying to query the 'SelectionStart' and
% 'SelectionEnd' properties.
%
% Basic options of entering user name and password with options for
% password validation.
%
% ToDo: Maybe add valid string input argument to only allow certain
% characters.

%% Program Information

% ProgramName = 'passwordEntryDialog';
% ProgramVersion = '0.0.2';
% svnRevision = '$Revision: 184 $';
% svnRevision = getSVNRevision(svnRevision);
% ProgramVersion = [ProgramVersion, '.' ,svnRevision];
%
% LastChangedDate = '$LastChangedDate: 2008-06-16 09:08:17 -0600 (Mon, 16 Jun 2008) $';
% ProgramDate = getSVNDate(LastChangedDate);

%% Check for Existance of Java
if ~usejava('swing')
   error('passwordEntryDialog: Java is required for this program to run.');
end

%% Parse Input Arguments
% Input arguments are parsed with the MATLAB inputParser class.

% Create input parser object
ProgramOptionsParser = inputParser;
ProgramOptionsParser.KeepUnmatched = true;

ProgramOptionsParser.addParamValue('enterUserName', false, @(x) islogical(x) || isnumeric(x));
ProgramOptionsParser.addParamValue('DefaultUserName', '', @ischar);
ProgramOptionsParser.addParamValue('ValidatePassword', false, @(x) islogical(x) || isnumeric(x));
ProgramOptionsParser.addParamValue('CheckPasswordLength', true, @(x) islogical(x) || isnumeric(x));
ProgramOptionsParser.addParamValue('PasswordLengthMin', 2, @isnumeric);
ProgramOptionsParser.addParamValue('PasswordLengthMax', 8, @isnumeric);
ProgramOptionsParser.addParamValue('WindowName', 'Login', @ischar);

% Parse Input Arguments
try
    ProgramOptionsParser.parse(varargin{:});
catch Error
    ProgramOptionsParser.parse;
    if strcmpi(Error.identifier, 'MATLAB:InputParser:ArgumentFailedValidation')
        error(Error.identifier, Error.message);
    end;
end;

ProgramOptions = ProgramOptionsParser.Results;

% Validate password length options
if ProgramOptions.CheckPasswordLength
    if ProgramOptions.PasswordLengthMax < ProgramOptions.PasswordLengthMin
        error('MATLAB:InputParser:ArgumentFailedValidation', 'PasswordLengthMax must be greater than PasswordLengthMin');
    end;
end;

%% Determine GUI Size and Position
% Center the GUI on the screen.

set(0,'Units','pixels')
Screen = get(0,'screensize');

PositionGUI = [0 0 300 50];
if ProgramOptions.enterUserName
    PositionGUI = PositionGUI + [0 0 0 50];
end;
if ProgramOptions.ValidatePassword
    PositionGUI = PositionGUI + [0 0 0 50];
    OffsetBottom = 43;
else
    OffsetBottom = 0;
end;

PositionGUI = [Screen(3:4)/2-PositionGUI(3:4)/2 PositionGUI(3:4)];
PositionLeft = 10;
BoxWidth = 200;

%% Create the GUI

BackgroundColor = 'w';
handles.figure1 = figure('Menubar','none', ...
    'Units','Pixels', ...
    'Resize','off', ...
    'NumberTitle','off', ...
    'Name',ProgramOptions.WindowName, ...
    'Position',PositionGUI, ...
    'Color', BackgroundColor, ...
    'WindowStyle','modal');

% Create Password Validation Entry Box
if ProgramOptions.ValidatePassword
    handles.java_PasswordValidate = javax.swing.JPasswordField();
    handles.java_PasswordValidate.setFocusable(true);
    [handles.java_PasswordValidate, handles.edit_PasswordValidate] = javacomponent(handles.java_PasswordValidate, [], handles.figure1);

    set(handles.edit_PasswordValidate, ...
        'Parent', handles.figure1, ...
        'Tag', 'edit_PasswordValidate', ...
        'Units', 'Pixels', ...
        'Position',[PositionLeft 10 BoxWidth 23]);

    handles.text_LabelPasswordValidate = uicontrol('Parent',handles.figure1, ...
        'Tag', 'text_LabelPassword', ...
        'Style','Text', ...
        'Units','Pixels',...
        'Position',[PositionLeft 33 BoxWidth 16], ...
        'FontSize',10, ...
        'String','Reenter password:',...
        'HorizontalAlignment', 'Left');
end;

% Create Password Entry Box
handles.java_Password = javax.swing.JPasswordField();
[handles.java_Password, handles.edit_Password] = javacomponent(handles.java_Password, [PositionLeft 10+OffsetBottom BoxWidth 23], handles.figure1);
handles.java_Password.setFocusable(true);

set(handles.edit_Password, ...
    'Parent', handles.figure1, ...
    'Tag', 'edit_Password', ...
    'Units', 'Pixels', ...
    'Position',[PositionLeft 10+OffsetBottom BoxWidth 23]);
drawnow;    % This drawnow is required to allow the focus to work
   
handles.text_LabelPassword = uicontrol('Parent',handles.figure1, ...
    'Tag', 'text_LabelPassword', ...
    'Style','Text', ...
    'Units','Pixels',...
    'Position',[PositionLeft 33+OffsetBottom BoxWidth 16], ...
    'FontSize',10, ...
    'String','Password:',...
    'HorizontalAlignment', 'Left');

% Create OK Pushbutton
handles.pushbutton_OK = uicontrol('Parent',handles.figure1, ...
    'Tag', 'pushbutton_OK', ...
    'Style','Pushbutton', ...
    'Units','Pixels',...
    'Position',[PositionLeft+BoxWidth+5 10 30 23], ...
    'FontSize',10, ...
    'String','OK',...
    'HorizontalAlignment', 'Center');

% Create Cancel Pushbutton
handles.pushbutton_Cancel = uicontrol('Parent',handles.figure1, ...
    'Tag', 'pushbutton_Cancel', ...
    'Style','Pushbutton', ...
    'Units','Pixels',...
    'Position',[PositionLeft+BoxWidth+30+7 10 50 23], ...
    'FontSize',10, ...
    'String','Cancel',...
    'HorizontalAlignment', 'Center');

% Create User Name Edit Box
if ProgramOptions.enterUserName
    handles.java_UserName = javax.swing.JTextField();
    handles.java_UserName.setFocusable(true);
    [handles.java_UserName, handles.edit_UserName] = javacomponent(handles.java_UserName, [], handles.figure1);

    set(handles.edit_UserName, ...
        'Parent', handles.figure1, ...
        'Tag', 'edit_UserName', ...
        'Units', 'Pixels', ...
        'Position',[PositionLeft 53+OffsetBottom 200 23]);
    set(handles.java_UserName, 'Text', ProgramOptions.DefaultUserName);
    drawnow;    % This drawnow is required to allow the focus to work

    handles.text_LabelUserName = uicontrol('Parent',handles.figure1, ...
        'Tag', 'text_LabelUserName', ...
        'Style','Text', ...
        'Units','Pixels',...
        'Position',[PositionLeft 76+OffsetBottom 200 16], ...
        'FontSize',10, ...
        'String','User name:',...
        'HorizontalAlignment', 'Left');

    %uicontrol(handles.edit_UserName);
    %set(handles.figure1,'CurrentObject',handles.java_UserName)
    handles.java_UserName.requestFocus;     % Get focus
    drawnow;
else
    handles.java_Password.requestFocus;     % Get focus
    drawnow;
end;

%% Setup Callbacks for Objects
% Adds the callback functions for the objects in the GUI

set(handles.pushbutton_OK,     'Callback', {@pushbutton_OK_Callback, handles, ProgramOptions}, 'KeyPressFcn', {@pushbutton_KeyPressFcn, handles, ProgramOptions});
set(handles.pushbutton_Cancel, 'Callback', {@pushbutton_Cancel_Callback, handles, ProgramOptions}, 'KeyPressFcn', {@pushbutton_KeyPressFcn, handles, ProgramOptions});
set(handles.java_Password, 'ActionPerformedCallback', {@pushbutton_OK_Callback, handles, ProgramOptions});

if ProgramOptions.ValidatePassword
    if ProgramOptions.enterUserName
        ObjectNext = handles.java_UserName;
    else
        ObjectNext = handles.java_Password;
    end;
    set(handles.java_PasswordValidate, 'ActionPerformedCallback', {@pushbutton_OK_Callback, handles, ProgramOptions}, 'NextFocusableComponent', ObjectNext);
    set(handles.java_Password, 'NextFocusableComponent', handles.java_PasswordValidate);
elseif ProgramOptions.enterUserName
    set(handles.java_Password, 'NextFocusableComponent', handles.java_UserName);
end;

if ProgramOptions.enterUserName
    set(handles.java_UserName, 'ActionPerformedCallback', {@pushbutton_OK_Callback, handles, ProgramOptions}, 'NextFocusableComponent', handles.java_Password);
end;

setappdata(handles.figure1, 'isCanceled', false);

%% Wait for Completion

% Wait for the user to complete entry.
drawnow;
uiwait(handles.figure1);

%% Get current information
% User either closed the window or pressed cancel or OK.

isCanceled = ~ishandle(handles.figure1) || getappdata(handles.figure1, 'isCanceled');
if isCanceled
    if ishandle(handles.figure1)
        delete(handles.figure1);
    end;
    Password = -1;
    UserName = '';
    return;
end;

Password = handles.java_Password.Password';
if ProgramOptions.enterUserName
    UserName = strtrim(get(handles.java_UserName, 'Text'));
else
    UserName = '';
end;
delete(handles.figure1);

%% DEFINE FUNCTIONS
% The subfunctions required by this program are in the following section.

function pushbutton_KeyPressFcn(hObject, eventdata, handles, ProgramOptions)

switch eventdata.Key
    case 'return'
        Callback = get(hObject, 'Callback');
        feval(Callback{1}, hObject, '', Callback{2:end});
end;

function pushbutton_OK_Callback(hObject, eventdata, handles, ProgramOptions)
if ProgramOptions.enterUserName
    % Check if username is blank
    UserName = strtrim(get(handles.java_UserName, 'Text'));
    if isempty(UserName)
        strMessage = 'Username is blank';
        %disp(strMessage)
        hError = errordlg(strMessage, 'passwordEntryDialog');
        uiwait(hError);
        return;
    end;
end;

if ProgramOptions.CheckPasswordLength
    %Password = handles.edit_Password.Password';
    Password = handles.java_Password.Password';
    if length(Password) < ProgramOptions.PasswordLengthMin || length(Password) > ProgramOptions.PasswordLengthMax
        strMessage = sprintf('Password must be between %d and %d characters', ...
            ProgramOptions.PasswordLengthMin, ...
            ProgramOptions.PasswordLengthMax);
        %disp(strMessage);
        hError = errordlg(strMessage, 'passwordEntryDialog');
        uiwait(hError);
        if ProgramOptions.ValidatePassword
            set(handles.java_PasswordValidate,'Text', '');
        end;
        handles.java_Password.requestFocus
        return;
    end;
end;

if ProgramOptions.ValidatePassword
    % Check if passwords match
    if ~isequal(handles.java_Password.Password, handles.java_PasswordValidate.Password)
        strMessage = 'Passwords do not match.  Please try again';
        %disp(strMessage);
        hError=errordlg(strMessage, 'passwordEntryDialog','modal');
        uiwait(hError);
        set(handles.java_Password,'Text', '');
        set(handles.java_PasswordValidate,'Text', '');

        handles.java_Password.requestFocus
        return;
    end;
end;
uiresume(handles.figure1);

function pushbutton_Cancel_Callback(hObject, eventdata, handles, ProgramOptions)
setappdata(handles.figure1, 'isCanceled', true);
uiresume(handles.figure1);
