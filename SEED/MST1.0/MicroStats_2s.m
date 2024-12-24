%MICROSTATS Calculates microstate statistics.
%
% Usage:
%   >> Mstats = MicroStats(X, A, L)
%   >> Mstats = MicroStats(X, A, L,polarity,fs)
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
%  Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
%  Inputs:
%   X - EEG (channels x samples (x trials)).
%   A - Spatial distribution of microstate prototypes (channels x  K).΢��״̬�ռ����
%   L - Label of the most active microstate at each timepoint (trials x
%       time).ÿ��ʱ��������Ծ��΢״̬�ı�ǩ
%
%  Optional input:
%   polarity - Account for polarity when fitting Typically off for
%              spontaneous EEG and on for ERP data (default = 0).
%   fs       - Sampling frequency of EEG (default = 1).
%
%  Outputs:
%  Mstats - Structure of microstate parameters per trial:ÿ�������΢��״̬�����ṹ
%   .Gfp        - Global field power
%   .Occurence  - Occurence of a microstate per s ÿ����ֵ�΢״̬
%   .Duration   - Average duration of a microstate ΢״̬��ƽ������ʱ��
%   .Coverage   - % of time occupied by a microstate ��΢��״̬ռ�õ�ʱ��
%   .GEV        - Global Explained Variance of microstate
%   .MspatCorr  - Spatial Correlation between template maps and microstates
%   .TP         - transition probabilities
%   .seq        - Sequence of reoccurrence of MS ((trials x)  time).
%   .msFirstTF  - First occurence of a microstate (similar to use like seq)
%                 ((trials x)  time)
%   .polarity   - see inputs.
%
%  Mstats.avgs - Structure of microstate parameters mean / and std over
%                trials.������΢��״̬������ֵ�ͱ�׼��Ľṹ��
%
% Authors:
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Z�rich, Psychologisches Institut, Methoden der
% Plastizit�tsforschung.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% September 2017.

% Copyright (C) 2017 Andreas Pedroni, andreas.pedroni@uzh.ch.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function Mstats = MicroStats_2s(X,A,L,polarity,fs)
%% Error check and initialisation ������ͳ�ʼ��
if nargin < 5
    fs = 1;
elseif nargin < 3
    polarity = 0;
elseif nargin < 3
    help MicroStats;
    return;
end

[C,N,Ntrials] = size(X); % Ҳ��˵���Ѷ�ά����������άΪ1����ά������Ҳ��ͬ���ǰ�nά����������n��1�ľ���һ��,��NtrialsΪ1
K = size(A,2);  % ���ص���A�������������΢״̬ԭ�͵���Ŀ


%% Spatial correlation between microstate prototypes and EEG ΢״̬ԭ�����Ե�ͼ�Ŀռ������
% force average reference
X = X - repmat(mean(X,1),[C,1,1]);  % mean(X,1)��ʾ���ط���X����ÿ�е�ƽ��ֵ1*samples��repmat(mean(X,1),[C,1,1])�����Ϊ62*samples

%Check if the data has more than one trial or not and reshape if necessary ��������Ƿ��в�ֹһ�ε�trail����Ҫʱ��������
if Ntrials > 1 % for epoched data
    X = squeeze(reshape(X, C, N*Ntrials));  % squeeze��������ɾ�������еĵ�һά���������һ��1x2x3�ľ���A��Ȼ��squeeze��A��������һ��2x3�ľ��󣬽���һάȴ��(��Ϊ��һλ��СΪ1)
end

% Normalise EEG and maps (average reference and gfp = 1 for EEG)
Xnrm = X ./ repmat(std(X,1), C, 1); % already have average reference �Ѿ�����ƽ���ο� std(X,1)����X�ı�׼�1��ʾ����N
A_nrm = (A - repmat(mean(A,1), C, 1)) ./ repmat(std(A,1), C, 1);

% Global map dissimilarity
GMD = nan(K,N);
for k = 1:K
    GMD(k,:) = sqrt(mean( (Xnrm - repmat(A_nrm(:,k),1,N)).^2 ));
end

% Account for polarity (recommended 0 for spontaneous EEG)
if polarity == 0
    GMDinvpol = nan(K,N);
    for k = 1:K
        GMDinvpol(k,:) = sqrt(mean( (Xnrm - repmat(-A_nrm(:,k),1,size(Xnrm,2))).^2));
    end
    idx = GMDinvpol < GMD;
    GMD(idx) = GMDinvpol(idx);
end

% Calculate the spatial correlation between microstate prototypes and EEG ����΢״̬ԭ�ͺ��Ե�ͼ֮��Ŀռ������
SpatCorr = 1 - (GMD.^2)./2;

% APedroni added this 13.12.2017
GEVtotal = sum((SpatCorr(sub2ind(size(SpatCorr),L,1:size(L,2))).*squeeze(std(X)).^2) ./sum(squeeze(std(X)).^2));

SpatCorr = squeeze(reshape(SpatCorr,K,N,Ntrials));

%% Sequentialize
% take into account the sequence of occurence of microstates. This makes
% only sense in ERP data (if it makes sense) ����΢��״̬�ĳ���˳����ֻ��ERP������������(���������Ļ�)
seq = nan(Ntrials, N);  % NtrialsΪ1��NΪsamples
msFirstTF = nan(Ntrials, N);
for trial = Ntrials
    first = 1;
    s = ones(K,1);  % ����5*1��ȫ1����
    for n = 1:(N-1)
        if L(trial,n) == L(trial,n+1)  % L��ʾlabel
            seq(trial,n) = s(L(trial,n));
            msFirstTF(trial,n) = first;
        else
            seq(trial,n) = s(L(trial,n))  ;
            s(L(trial,n),1) = s(L(trial,n),1) + 1;
            first = n;
            msFirstTF(trial,n) = first;
        end
    end
    seq(trial,n+1) = seq(trial,n);
    msFirstTF(trial,n+1) = msFirstTF(trial,n);
end

%% ��������trail��MspatCorr��GEV
MspatCorr = nan(Ntrials,K);
GEV = nan(Ntrials,K);
GFP = squeeze(std(X));  % std(X)��������������ı�׼��,��ʱ���Ե���N-1��Ĭ��Ϊ���У�squeezeȥ��sizeΪ1��ά��
if Ntrials>1
    GFP = GFP';
end
for k = 1:K
    for trial = 1:Ntrials
         % Average spatial correlation per Microstate ÿ��΢״̬��ƽ���ռ������
        MspatCorrTMP = SpatCorr(:,:,trial);
        MspatCorr(trial,k) = nanmean(MspatCorrTMP(k,L(trial,:)==k));
        
        % global explained variance Changed by Pedroni 3.1.2018
        GEV(trial,k) = sum( (GFP(trial,L(trial,:)==k) .* MspatCorrTMP(k,L(trial,:)==k)).^2) ./ sum(GFP(trial,:).^2);
    end
end

%% ������trail����������2s�����ݼ���4��ʱ��ͳ�Ʋ���
list = size(L,2);
if list == 47001
    Lcut = L(:,201:47000);
    Xcut = X(:,201:47000);
elseif list == 46601
    Lcut = L(:,201:46600);
    Xcut = X(:,201:46600);
elseif list == 41201
    Lcut = L(:,1:41200);
    Xcut = X(:,1:41200);
elseif list == 47601
    Lcut = L(:,1:47600);
    Xcut = X(:,1:47600);
elseif list == 37001
    Lcut = L(:,201:37000);
    Xcut = X(:,201:37000);
elseif list == 39001
    Lcut = L(:,201:39000);
    Xcut = X(:,201:39000);
elseif list == 47401
    Lcut = L(:,201:47400);
    Xcut = X(:,201:47400);
elseif list == 43201
    Lcut = L(:,1:43200);
    Xcut = X(:,1:43200);
elseif list ==53001
    Lcut = L(:,201:53000);
    Xcut = X(:,201:53000);
end
m = size(Lcut,2);
s = 1;
for i = 1:400:m
    sampleL = Lcut(:,i:i+399);
    sampleX = Xcut(:,i:i+399);
%% Microstate order for transition probabilities ΢״̬���е�ԾǨ����
order = cell(Ntrials,1);
for trial = 1:Ntrials
    [order{trial}, ~ ] = my_RLE(sampleL(trial,:));
end


%% Preallocating arrays for stats and readying GFP
% prepare arrays (only MGFP can have NANs!)
MGFP = nan(Ntrials,K);  % ����1*5�ľ���
MDur = zeros(Ntrials,K);
MOcc = zeros(Ntrials,K);
TCov = zeros(Ntrials,K);

GFPcut = squeeze(std(sampleX));  % std(X)��������������ı�׼��,��ʱ���Ե���N-1��Ĭ��Ϊ����
if Ntrials>1
    GFPcut = GFPcut';
end

%% For each MS class...
for k = 1:K
    for trial = 1:Ntrials
        
        % Mean GFP
        MGFP(trial,k) = nanmean(GFPcut(trial,sampleL(trial,:)==k));  % �����Ծ���GFPÿһ�����ֵ������֮ǰ����A�е�NaNԪ��
        
        [runvalue, runs] = my_RLE(sampleL(trial,:));
        
        % Mean Duration
        if isnan(nanmean(runs(runvalue == k)))
            MDur(trial,k) = 0;
            MOcc(trial,k) = 0;
        else
            MDur(trial,k) =  nanmean(runs(runvalue == k)) .* (1000 / fs);
            % Occurence
            MOcc(trial,k) =  length(runs(runvalue == k))./400.* fs;
        end
        % time coverage
        TCov(trial,k) = (MDur(trial,k) .* MOcc(trial,k))./ 1000;
        
       
    end
end
MGFPsample(s,:) = MGFP;
MDursample(s,:) = MDur;
MOccsample(s,:) = MOcc;
TCovsample(s,:) = TCov;
%% Transition Probabilities (as with hmmestimate(states,states);)
for trial = 1:Ntrials
    states = order{trial,:};
    states = states(states ~= 0);
    % prepare output matrix
    numStates = K;
    tr = zeros(numStates);
    seqLen = length(states);
    % count up the transitions from the state path
    for count = 1:seqLen-1
        tr(states(count),states(count+1)) = tr(states(count),states(count+1)) + 1;
    end
    trRowSum = sum(tr,2);
    % if we don't have any values then report zeros instead of NaNs.
    trRowSum(trRowSum == 0) = -inf;
    % normalize to give frequency estimate.
    TP(:,:,trial) = tr./repmat(trRowSum,1,numStates);
end
TPsample{s,:} = TP;
s = s+1;
end

%% Write to EEG structure
%   per trial:
Mstats.GEVtotal = GEVtotal;
Mstats.Gfp = MGFPsample;
Mstats.Occurence = MOccsample;
Mstats.Duration = MDursample;
Mstats.Coverage = TCovsample;
Mstats.GEV = GEV;
Mstats.MspatCorr = MspatCorr;
Mstats.TP = TPsample;
Mstats.seq = seq;
Mstats.msFirstTF = msFirstTF;
Mstats.polarity = polarity;

if Ntrials > 1
    % mean parameters
    Mstats.avgs.Gfp = nanmean(MGFP,1);
    Mstats.avgs.Occurence = nanmean(MOcc,1);
    Mstats.avgs.Duration = nanmean(MDur,1);
    Mstats.avgs.Coverage = nanmean(TCov,1);
    Mstats.avgs.GEV = nanmean(GEV,1);
    Mstats.avgs.MspatCorr = nanmean(MspatCorr,1);
    
    % standard deviation of parameters ������׼��
    Mstats.avgs.stdGfp = nanstd(MGFP);
    Mstats.avgs.stdOccurence = nanstd(MOcc);
    Mstats.avgs.stdDuration = nanstd(MDur);
    Mstats.avgs.stdCoverage = nanstd(TCov);
    Mstats.avgs.stdGEV = nanstd(GEV);
    Mstats.avgs.stdMspatCorr = nanstd(MspatCorr);
end
end

function [d,c]=my_RLE(x)
%% RLE
% This function performs Run Length Encoding to a strem of data x.�ú�����������xִ�����г��ȱ���
% [d,c]=rl_enc(x) returns the element values in d and their number of  ����d�е�Ԫ��ֵ������ִ�������c�У�
% apperance in c. All number formats are accepted for the elements of x.  x�е�Ԫ�ؿ��Խ����������ָ�ʽ��
% This function is built by Abdulrahman Ikram Siddiq in Oct-1st-2011 5:15pm.

if nargin~=1
    error('A single 1-D stream must be used as an input') % ����ʹ�õ���һά����Ϊ����
end

ind=1;
d(ind)=x(1);
c(ind)=1;

for i=2 :length(x)
    if x(i-1)==x(i)
        c(ind)=c(ind)+1;
    else ind=ind+1;
        d(ind)=x(i);
        c(ind)=1;
    end
end

end
