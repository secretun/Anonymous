% 获取E:\DEAP\mat目录下所有.mat文件
files = dir('E:\DEAP\mat\*.mat');
% 遍历每个文件
for i = 1:length(files)
    % 加载.mat文件
    filePath = fullfile(files(i).folder, files(i).name);
    dataStruct = load(filePath);
    % 假设data字段存在，并且是40x40x8064的数组
    % 提取前32列，得到40x32x8064的数组
    if isfield(dataStruct, 'data') && size(dataStruct.data, 2) == 40
        extractedData = dataStruct.data(:, 1:32, :);
        % 定义保存路径
        savePathData = fullfile('E:\DEAP\mat_data\', files(i).name);
        % 保存提取后的数据
        save(savePathData, 'extractedData');
    else
        disp(['Data field is missing or does not have the expected dimensions in ', files(i).name]);
    end
end

