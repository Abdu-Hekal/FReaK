function setupKF()

try
    %go to koopman falsification home folder
    cdr = pwd; %current folder
    filePath = which('setupKF'); %filepath
    [kfFolder, ~, ~] = fileparts(filePath);
    cd(kfFolder);
    % add to path and remove archive and hacky cora files
    addpath(genpath(kfFolder))
    rmpath(fullfile('external','CORA'))
    rmpath(genpath('archive'))

    %check installed toolbox
    %code taken from CORA installation.
    auxInstallToolbox('Symbolic Math Toolbox');

    %install mpt and yalmip
    installMPT()

    %install python libraries
    pythonLibInstall();

    %make folder to store auxilary tools
    folderName = 'auxilary';
    if ~(exist(folderName, 'dir') == 7)
        mkdir(folderName);
    end

    %install breach and CORA
    if (exist('InitBreach','file')==2)
        disp('Breach already installed.');
    else
        userResponse = input('Breach not found, is it installed? (y/n)', 's');
        if userResponse
            error('please add Breach base folder with (InitBreach.m) to path')
        else
            zipURL = 'https://github.com/decyphir/breach/archive/refs/heads/master.zip';
            installExtToolbox('breach',zipURL,'breach.zip')
            installBreach; %compiles C/C++ mex functions used by Breach
        end
    end
    if (exist('test_requiredToolboxes','file')==2)
        disp('CORA already installed.');
    else
        userResponse = input('CORA not found, is it installed? (y/n)', 's');
        if userResponse
            error('please add CORA base folder with (test_requiredToolboxes.m) to path')
        else
            zipURL = 'https://tumcps.github.io/CORA/data/archive/version/CORA_2022.zip';
            installExtToolbox('CORA_2022',zipURL,'CORA_2022.zip');
        end
    end

    %setup Breach and CORA
    InitBreach()
    setupBreach()
    setupCora()

    %setup gurobi
    setupGurobi()

    %go back to base folder
    cd(cdr);

    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    disp('Koopman Falsification Successfully Setup!')
    disp('Run "savepath;" to store path setup for future matlab sessions')

catch ME
    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    fprintf(2,'Setup Failed due to following reason, please fix and rurun setup: \n')
    disp(ME.message)
end

end

% -------------------------- Auxiliary Functions --------------------------
function res = auxInstallToolbox(text)
res = auxTestToolboxInstallation(text);
while ~res
    auxDisplayInstallPrompt(text);
    res = auxTestToolboxInstallation(text);
end
end
function res = auxTestToolboxInstallation(text)
res = any(any(contains(struct2cell(ver),text)));
end
function auxDisplayInstallPrompt(text)
error(['''<strong>%s</strong>'' is missing and requires manual installation. \n' ...
    '  Please install it via the MATLAB Add-On Explorer. \n'], text)
end

function pythonLibInstall()
%ensure python is installed
assert(~isempty(pyenv), 'Python is not installed or configured in MATLAB.')
%check if correct python environment is setup
disp(pyenv)
userResponse='';
while ~strcmpi(userResponse, 'y') && ~strcmpi(userResponse, 'n')
    userResponse = input('Is this the correct python environment? (y/n): ', 's');
end
if strcmpi(userResponse, 'n')
    if pyenv().Status == "Loaded"
        error('To change the Python version, restart MATLAB and run setupKF again')
    end
    userResponse = input('please enter executable of python environment to be used: ', 's');
    pyenv(Version=userResponse)
end

%get pip for python executable
pythonExecutablePath = char(pyenv().Executable);
% Split the path into components
[baseFolderPath, ~, ~] = fileparts(pythonExecutablePath);
% Specify the new last part of the path
newLastPart = 'pip';
% Construct the new path
pipExecutablePath = fullfile(baseFolderPath, newLastPart);

%install autokoopman
% Check if autokoopman is installed
try
    py.importlib.import_module('autokoopman');
    disp('Autokoopman already installed.')
catch
    [status,result]=system([pipExecutablePath ' install autokoopman']);
    assert(status == 0, ['Installation of autokoopman failed. Error message: ', result]);
    disp('Successfully installed Autokoopman!')
end
end

function installExtToolbox(name,zipURL,zipName)
% Specify the local folder to save the ZIP file
destinationFolder = 'auxilary';  % Replace with the desired local folder name
% Check if the destination folder already exists
if exist(fullfile(destinationFolder, name), 'dir') == 7
    disp([name ' already downloaded.']);
else
    % Download the ZIP file
    try
        zipFileName = fullfile(destinationFolder, zipName);
        webopt=weboptions('Timeout',30); %timeout of 30s
        websave(zipFileName, zipURL,webopt);
        % Unzip the downloaded file
        unzip(zipFileName, fullfile(destinationFolder,name));
        %add to path
        addpath(genpath(fullfile(destinationFolder,name)))
        % Delete the ZIP file
        delete(zipFileName);

        disp(['Succesfully downloaded ' name]);
    catch ME
        disp(['Failed to download ' name ', please manually download and add to path'])
        rethrow(ME)
    end

end
end

function setupBreach()
%remove conflicting breach files from path
filePath = which('InitBreach');
[breachFolder, ~, ~] = fileparts(filePath);
WarnState=warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath(fullfile(breachFolder, 'Ext')))
rmpath(genpath(fullfile(breachFolder, 'Examples')))
rmpath(genpath(fullfile(breachFolder, 'Online')))
warning(WarnState);
end

function setupCora()
filePath = which('stl'); %cora stl filepath
[coraStlFolder,~,~] = fileparts(filePath);

%replace CORA stl files to add string representation and abs function, and
%fix negationNormForm
sourceFolder = fullfile('external','CORA');
destinationFolder = coraStlFolder;
files = '*.m';  % all m files
% Construct the full paths for the source and destination files
sourceFiles = fullfile(sourceFolder, files);
% Move the files with the option to overwrite existing files
[success,message]=copyfile(sourceFiles, destinationFolder, 'f');
assert(success, ['File copy failed with message ' message]);
%remove CORA modification files from path
if contains(path, sourceFolder)
    rmpath(sourceFolder)
end
end

function setupGurobi()
userResponse='';
while ~strcmpi(userResponse, 'y') && ~strcmpi(userResponse, 'n')
    userResponse = input('Koopman falsification depends on optimization solver. Examples use Gurobi,  \n do you want to check for gurobi setup? (y/n): ', 's');
end
if strcmpi(userResponse, 'y')
    status = system('gurobi_cl');
    if status==0
        disp('Valid Gurobi license found!')
    else
        error('No valid Gurobi license found')
    end
    if exist('gurobi_setup','file')
        gurobi_setup;
        disp('')
        disp('Gurobi Succesfully setup!')
    else
        error('Gurobi for matlab is not downloaded, follow instructions here: <a href="https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab">Gurobi Installation Guide for MATLAB</a>')
    end
end
end


