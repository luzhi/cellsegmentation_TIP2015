clear;
close all;
clc;
warning off;
pause(0.5);

% **********************************************
% Need to load STPRTOOL & VLfeat libraries here!
% **********************************************
% STPRTOOL: http://cmp.felk.cvut.cz/cmp/software/stprtool/
% VLfeat: http://www.vlfeat.org/

%==================================
%   Add folders/subfolders under
%   current folder
%==================================
currentfolder = pwd;
addpath(genpath(currentfolder))

%+------------------------+
%| Distance Map Parameter |
%+------------------------+
beta_logistic_set = (5)';

%+---------------------------+
%| Joint Level Set Parameter |
%+---------------------------+
kappa_set = (13)';

chi_set = (3)';

%+---------------------------------+
%| Joint Level Set Updating Policy |
%+---------------------------------+
loop = 5;

storageCommonPath = {'Common/'};
storageInitial= {'Initial/'};
storageDist= {'LSF/'};

DSOption = 'EDF';
%DSOption = 'SyntheticTrain';
%DSOption = 'SyntheticTest';

%+----------------+
%| Create folders |
%+----------------+
if strcmp(DSOption, 'EDF')
    OUTPUT_PATH = {'OUT_EDF/'};
end

if strcmp(DSOption, 'SyntheticTrain')
    OUTPUT_PATH = {'OUT_SynTrain/'};
end

if strcmp(DSOption, 'SyntheticTest')
    OUTPUT_PATH = {'OUT_SynTest/'};
end


mkdir(OUTPUT_PATH{:});
mkdir(strcat(OUTPUT_PATH{:}, storageCommonPath{:}));
mkdir(strcat(OUTPUT_PATH{:}, storageInitial{:}));
mkdir(strcat(OUTPUT_PATH{:}, storageDist{:}));
for i = 1:loop
    mkdir(strcat(OUTPUT_PATH{:}, storageDist{:}, 'LSF', num2str(i), '/'));
end

%+-------------------------------------+
%| Run Algorithm for Different DataSet |
%+-------------------------------------+
if strcmp(DSOption, 'EDF')
    %+---------------------------+
    %|      Pap smear - EDF      |
    %+---------------------------+
    Runner_inOne(DSOption, './ims/EDF/','./ims/EDF/', ...
                 storageCommonPath{:}, storageInitial{:}, storageDist{:}, ...
                 beta_logistic_set, kappa_set, chi_set, loop, ...
                 OUTPUT_PATH{:})
end

if strcmp(DSOption, 'SyntheticTest') || strcmp(DSOption, 'SyntheticTrain')
    %+-------------------------------+
    %| Pap smear - New Train Dataset |
    %+-------------------------------+
    Runner_inOne(DSOption, './ims/Synthetic/','./ims/Synthetic/', ...
                 storageCommonPath{:}, storageInitial{:}, storageDist{:}, ...
                 beta_logistic_set, kappa_set, chi_set, loop, ...
                 OUTPUT_PATH{:})
end
