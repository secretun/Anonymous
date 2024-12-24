% 定义输入输出路径
setDir = 'E:\DEAP\sub\sub9\分好trail无通道';
set1Dir = 'E:\DEAP\sub\sub9\分好trail有通道';
chanFilePath = 'E:\DEAP\chan\channel_32.mat';
% 检查 set1 目录是否存在，如果不存在则创建
if ~exist(set1Dir, 'dir')
    mkdir(set1Dir);
end
% 加载电极信息
chanData = load(chanFilePath);
chanlocs = chanData.EEG.chanlocs;
% 获取 set 目录下的所有 .set 文件
setFiles = dir(fullfile(setDir, '*.set'));
eeglab;
% 遍历所有 .set 文件
for i = 1:length(setFiles)
    % 加载 .set 文件
    setFilePath = fullfile(setDir, setFiles(i).name);
    EEG = pop_loadset('filename', setFilePath);
    % 添加电极信息
    EEG.chanlocs = chanlocs;
    EEG = eeg_checkset(EEG);
    % 定义新的文件名和保存路径
    [~, name, ~] = fileparts(setFiles(i).name);
    newSetFilePath = fullfile(set1Dir, [name, '.set']);
    % 保存更改到新目录
    EEG = pop_saveset(EEG, 'filename', newSetFilePath);
end