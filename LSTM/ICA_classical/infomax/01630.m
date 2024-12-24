load d:\important\eegdat\eeg1.mat
wc=s1base1;
[w,s]=runica(wc);
act=icaact(wc,w,s,0);
[pro] = plotproj(wc,w,s,1);
 [chandata] = chanproj(projdata,1,10);
fprintf('\n Note that two components make up most of the response peak.\n');
