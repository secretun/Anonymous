clc
clear all
%% 读取数据整合
[feature_data, txt, raw2] = xlsread('E:\DEAP\sub\cross_fea.xlsx',6,'A1:GS2394');
data = feature_data(:,2:201);
data(isnan(data)) = 0;
[data_sta, ps] = mapminmax(data',-1,1);%数据归一化
data = data_sta';
len = size(data,1);
label = zeros(len,1);
label(1:978) = 1;
label(979:1998) = 2;
label(1999:1394) = 3;
data = data(:,any(data));%去除全为0的列
%% 整体数据分成训练集和测试集
indices = crossvalind('Kfold', max(size(data,1)), 10);%将数据样本随机分割为10部分，将数据集划分为10个大小相同的互斥子集
svmpredictlable = cell(10,1);%SVM
knnpredictlabel = cell(10,1);%KNN
accuracy_svm = zeros(10,3);
accuracy_knn = zeros(10,1);
for i = 1:1:10
test        = (indices == i);
train       = ~test;%1：表示该组数据被选中，0：未被选中
traindata   = data(train, :);
testdata    = data(test, :);
train_label = label(train,:);%label数组存放情感的三种分类情况
test_label  =label(test,:);
 %%随机排列
newtrainX=[];newtrainY=[];newtestX=[];newtestY=[];
perm1=randperm(length(traindata(:,1)));
newtrainX(:,:)=traindata(perm1,:);
newtrainY(:,:)=train_label(perm1,:);
perm2=randperm(length(testdata(:,1)));
newtestX(:,:)=testdata(perm2,:);
newtestY(:,:)=test_label(perm2,:);
%% SVM
 c = 10;
 g = 0.5;
model = fitcecoc(newtrainX, newtrainY, 'Learners', templateSVM('KernelFunction', 'rbf', 'BoxConstraint', c, 'KernelScale', 1/g));
[svmpredict_label, ~] = predict(model, newtestX);
accuracys = sum(svmpredict_label == newtestY) / length(newtestY) * 100;
accuracy_svm(i,:) = accuracys;
svmpredictlable{i,1} = svmpredict_label;
    %% KNN测试
% [knnpredict_label] = KNN(newtrainX,newtrainY,newtestX, 6);
% numCorrect = sum(test_label == knnpredict_label);  % 统计正确分类的样本数量
% numTotal = numel(test_label);  % 总样本数量
% accuracyk = (numCorrect / numTotal) * 100;  % 计算百分比准确率
% %[corrPredictions, accuracyk] =  Misclassification_accuracy(newtestY, knnpredict_label);
% accuracy_knn(i,:) = accuracyk;
% knnpredictlabel{i,1} = knnpredict_label;
end
mean_accuracys =mean(accuracy_svm(:,1));
svm_sd=std(accuracy_svm(:,1));
%  mean_accuracyk = mean(accuracy_knn);
%  knn_sd=std(accuracy_knn);