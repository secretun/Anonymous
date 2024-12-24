% 定义输入输出路径
matDir = 'E:\DEAP\mat_data';
setDir = 'E:\DEAP\set';

% 检查输出目录是否存在，如果不存在则创建
if ~exist(setDir, 'dir')
    mkdir(setDir);
end

% 获取 mat 目录下的所有 .mat 文件
matFiles = dir(fullfile(matDir, '*.mat'));

% 初始化EEGLAB，确保EEGLAB函数可用
eeglab;

% 遍历所有 .mat 文件
for i = 1:length(matFiles)
    % 加载 .mat 文件
    matFilePath = fullfile(matDir, matFiles(i).name);
    data = load(matFilePath, 'extractedData'); % 加载extractedData字段

    % 创建EEG结构体
    EEG = eeg_emptyset();
    EEG.data = data.extractedData;
    EEG.srate = 128; % 设置采样率
    EEG = eeg_checkset(EEG);

    % 应用平均参考
    EEG = pop_reref(EEG, []);

    % 定义新的文件名和保存路径
    [~, name, ~] = fileparts(matFiles(i).name);
    setFilePath = fullfile(setDir, [name, '.set']);

    % 保存为 .set 文件
    EEG = pop_saveset(EEG, 'filename', setFilePath);
end


