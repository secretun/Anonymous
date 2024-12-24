clc;clear;
% 定义源文件和目标文件的路径
sourceFile = 'E:\SEED\microstate\分开导入\实验融合\9_15.xlsx';
targetFile = 'E:\SEED\microstate\分开导入\实验融合\fea1.xlsx';
% 循环处理每个工作表
for sheet = 1:6
    % 读取源文件中的数据
    sourceData = readtable(sourceFile, 'Sheet', sheet);
    % 如果目标文件的工作表已经存在数据，则读取这些数据
    if exist(targetFile, 'file') == 2
        try
            targetData = readtable(targetFile, 'Sheet', sheet);
        catch
            targetData = [];
        end
    else
        targetData = [];
    end
    % 将源文件中的数据追加到目标文件的数据后面
    newData = [targetData; sourceData];
    % 将新的数据写入目标文件的工作表，不写入变量名
    writetable(newData, targetFile, 'Sheet', sheet, 'WriteVariableNames', false);
end

