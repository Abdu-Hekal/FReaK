function autoSetupKF()

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

%setup Breach and CORA
InitBreach()
setupBreach()
setupCora()

%setup mpt
mpt_init

%add gurobi to path
addpath(genpath('/opt/gurobi/linux64/matlab'));
% run setup
gurobi_setup

%go back to base folder
cd(cdr);
%save modified path
savepath;

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp('Koopman Falsification Successfully Setup!')

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

function setupBreach()
%remove conflicting breach files from path
filePath = which('InitBreach');
[breachFolder, ~, ~] = fileparts(filePath);
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath(fullfile(breachFolder, 'Ext')))
rmpath(genpath(fullfile(breachFolder, 'Examples')))
warning('on', 'MATLAB:rmpath:DirNotFound');
end

function setupCora()
filePath = which('stl'); %cora stl filepath
[coraStlFolder,~,~] = fileparts(filePath);

%replace CORA stl files to add string representation and abs function
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


