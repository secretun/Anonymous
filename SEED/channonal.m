% 设置文件夹路径
inputFolderPath = 'E:\SEED\SEED_EEG\nochanol\15_20131105';
outputFolderPath = 'E:\SEED\SEED_EEG\yeschanol\15_20131105';
% 获取文件夹中的所有 .set 文件
setFiles = dir(fullfile(inputFolderPath, '*.set'));
% 初始化 EEGLAB 变量
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% 遍历文件并应用操作
for i = 1:numel(setFiles)
    % 构建输入文件的完整路径
    inputFilePath = fullfile(inputFolderPath, setFiles(i).name);
    % 加载 EEG 数据集
    EEG = pop_loadset('filename', setFiles(i).name, 'filepath', inputFolderPath);
    % 编辑通道位置信息
%     EEG = pop_chanedit(EEG, 'lookup', 'E:\SEED\SEED_EEG\官网SEED\channel_62_pos.locs');
%     EEG=pop_chanedit(EEG, 'load',{'E:\SEED\SEED_EEG\官网SEED\channel_62_pos.locs','filetype','autodetect'});  
    EEG.srate=128;
    EEG = pop_reref(EEG, []); %设置平均参考
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
    % 构建输出文件的完整路径
    outputFilePath = fullfile(outputFolderPath, setFiles(i).name);
    % 保存处理后的 EEG 数据集到指定文件夹
    EEG = pop_saveset(EEG, 'filename', setFiles(i).name, 'filepath', outputFolderPath);
   clear EEG
end
% 在完成所有文件的处理后，重绘 EEGLAB 界面
eeglab redraw;
