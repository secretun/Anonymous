% 定义.mat文件的源目录
sourceDir = 'E:\DEAP\mat_data\三维';
% 定义输出目录为二维数组的保存位置
outputDir = 'E:\DEAP\mat_data\二维';
% 获取sourceDir目录下所有.mat文件
files = dir(fullfile(sourceDir, '*.mat'));
% 遍历每个文件
for i = 1:length(files)
    % 加载.mat文件
    filePath = fullfile(files(i).folder, files(i).name);
    dataStruct = load(filePath);
    % 假设数据存储在名为data的字段中
    if isfield(dataStruct, 'data')
        data = dataStruct.data; % 提取三维数组数据 
        % 将三维数组转换为二维数组
        [dim1, dim2, dim3] = size(data);
        data2D = reshape(data, dim1, dim2*dim3);
        % 保存转换后的二维数组，保持文件名不变
        savePath = fullfile(outputDir, files(i).name);
        save(savePath, 'data2D');
    else
        warning('No data field found in %s.', files(i).name);
    end
end
