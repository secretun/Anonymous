% 定义源目录和目标目录
sourceDirs = {'E:\DEAP\lable\1', 'E:\DEAP\lable\2'};
targetDirs = {'E:\DEAP\lable\1\1_1', 'E:\DEAP\lable\2\2_1'};
labels = {'label1', 'label2'};
% 遍历每个源目录
for i = 1:length(sourceDirs)
    % 获取当前目录下所有.mat文件
    files = dir(fullfile(sourceDirs{i}, '*.mat'));
    
    % 确保目标目录存在
    if ~exist(targetDirs{i}, 'dir')
        mkdir(targetDirs{i});
    end
    
    % 遍历每个文件
    for j = 1:length(files)
        % 加载.mat文件
        filePath = fullfile(files(j).folder, files(j).name);
        data = load(filePath);
        
        % 检查并修改对应的label字段
        if isfield(data, labels{i})
            labelData = data.(labels{i});
            fprintf('Original data in %s: %s\n', filePath, mat2str(labelData));
            
            % 修改数值：1~3的为-1，4~6为0，7~9为1
            labelData(labelData >= 1 & labelData <= 3) = -1;
            labelData(labelData >3 & labelData <= 6) = 0;
            labelData(labelData >6 & labelData <= 9) = 1;

            data.(labels{i}) = labelData;
            
            fprintf('Modified data: %s\n', mat2str(labelData));
            
            % 保存修改后的数据到新的目录
            savePath = fullfile(targetDirs{i}, files(j).name);
            save(savePath, '-struct', 'data');
        else
            warning('Field %s not found in %s.', labels{i}, files(j).name);
        end
    end
end
