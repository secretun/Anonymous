function [result,new_label,num_seg] = SEEDSingleSubData1(sub_no,chan_select,seg_len,overlap,ica,max_step,lrate,add_noise,noise_range)

addpath I:\SEED\Preprocessed_EEG
fileFolder = fullfile('I:\SEED\Preprocessed_EEG');
dirOutput = dir(fullfile(fileFolder,'*mat'));
allfilenames = {dirOutput.name};
subfilename = allfilenames(sub_no);
% subfilename = allfilenames(3*sub_no-2:3*sub_no);
num_seg = [];
result = [];

if chan_select
    chan_idx = [1,3,4,5,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,53,55,59,60,61];
%     chan_idx = [1,2,3,15,23,24,32,33,41]; 
%     chan_idx = [15,23,24,32,33,41];
else
    chan_idx = [1:1:62];
end

l = length(subfilename);
for i = 1:l
    filename = subfilename{i};
    data = load(filename);
    fields = fieldnames(data);
    for j = 1:length(fields)
        X = getfield(data,fields{j});
        X = X(chan_idx,:);
        [result1,num_seg1] = DataSeg(X,seg_len,overlap,ica,max_step,lrate,add_noise,noise_range);
        result = cat(3,result,result1);
        num_seg = cat(2,num_seg,num_seg1);
    end
end

label = load('I:\SEED\label.mat');
label = label.label;
label = repmat(label,1,l);
new_label = [];
for i = 1:length(label)
    label0 = label(i) * ones(1,num_seg(i));
    new_label = [new_label,label0];
end
end
