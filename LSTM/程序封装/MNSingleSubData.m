function [result,new_label] = MNSingleSubData(sub_no,seg_len,overlap,ica,max_step,lrate,add_noise,noise_range)

filenames_mahnob = load('valence_filenames.mat');
filenames_mahnob = filenames_mahnob.filename;
nt_mahnob = length(filenames_mahnob);

for l1 = 1:nt_mahnob
    filenames_mahnob{1,l1}(1) = 'G';
    filenames_mahnob{1,l1}(10:22) = '';
end

file_sub = ['Part_',num2str(sub_no),'_S_Trial'];
file_index = find(contains(filenames_mahnob,file_sub)==1);
num_seg = [];
result = [];
for i = file_index
    file = filenames_mahnob{i};
    [~,X] = edfread(file);
    X = X(1:32,:);
    [result1,num_seg1] = DataSeg(X,seg_len,overlap,ica,max_step,lrate,add_noise,noise_range);
    result = cat(3,result,result1);
    num_seg = cat(2,num_seg,num_seg1);
end

label = load('valence_label.mat');
label = label.valence_label;
label = label(file_index);
new_label = [];
for i = 1:length(label)
    label0 = label(i) * ones(num_seg(i),1);
    new_label = [new_label;label0];
end

end