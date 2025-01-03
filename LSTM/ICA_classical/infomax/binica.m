% binica - Run stand-alone binary version of runica() from the
%          Matlab command line. Saves time and memory. 
%          If stored in a file, data are not read into Matlab.
% Usage:
%         >> [wts,sph] = binica(datavar, [ica keywords, args]);
%  Else   >> [wts,sph] = binica('datafile',chans,frames,[ica_keys,args]);
%
% Inputs:
%         datavar = (chans,times) data matrix 
%         datafile = quoted 'filename' of float data file
%
% Allowed ica_keys and default args:
%
%         extended       0   [default: don't look for subGaussian components]
%         pca            0   [default: don't reduce data dimension]
%         blocksize      0   [default: heuristic dependent on data size]
%
%         lrate          1e-4    %        stop           1e-6
%         maxsteps       512     %        annealstep     0.98  [range 0-1]
%         annealdeg      60      %        momentum       0     [range 0-1]
%         sphering       on      %        bias           on
%         posact         on      %        verbose        on

% Scott Makeig 8/4/00 - CNL / The Salk Institute, La Jolla CA
% Calls ica implementation of Sigurd Enghoff, translation of runica()

% 08/07/00 Added warning to update icadefs.m -sm
% 09/08/00 Added tmpint to script, weights and sphere files to avoid
%          conflicts when multiple binica sessions run in the same pwd -sm
% 10/07/00 Fixed bug in reading wts when 'pca' ncomps < nchans -sm
 
function [wts,sph] = binica(data,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,var21,var22,var23,var24,var25)

if nargin < 1 | nargin > 25
    more on
    help binica
    more off
    return
end

icadefs % import ICA and SC
if ~exist('SC')
  fprintf('binica: You need to update your icadefs file to include var SC.\n')
  return
end

%
% select random integer 1-10000 to index the binica data files
% make sure no such script file already exists in the pwd
%
scriptfile = SC;
while exist(scriptfile)
  tmpint = randperm(10000);
  tmpint = int2str(tmpint(10000));
  scriptfile = ['binica' tmpint '.sc'];
end

nchans = 0;
tmpdata = [];
if ~isstr(data) % data variable name given
  if ~exist('data')
    fprintf('\nbinica(): Variable name data not found.\n');
    return
  end
  nchans = size(data,1);
  nframes = size(data,2);
  tmpdata = ['binica' tmpint '.fdt'];
  floatwrite(data,tmpdata);
  datafile = tmpdata;
  firstarg = 2;
else % data filename given
  if ~exist(data)
    fprintf('\nbinica(): File data not found.\n')
    return
  end
  datafile = data;
  if nargin<3
    fprintf(...
'\nbinica(): Data file name must be followed by chans, frames\n');
    return
  end
  nchans = var2;
  nframes = var3;
  if isstr(nchans) | isstr(nframes)
    fprintf(...
'\nbinica(): chans, frames args must be given after data file name\n');
    return
  end
  firstarg = 4;
end
%
% read in the master ica script file SC
%
flags = [];
args  = [];
fid = fopen(SC,'r');
if fid < 0
  fprintf('\nbinica(): Master .sc file %s not read!\n',SC)
     return
end

%
% read SC file info into flags and args lists
%
s = [];
f = 0; % flag number in file
while isempty(s) | s ~= -1
 s = fgetl(fid);
  if s ~= -1
   if ~isempty(s)
     s = rmcomment(s,'#');
     if ~isempty(s)
       f = f+1;
       s = rmspace(s);
       [w s]=firstword(s);
       if ~isempty(s)
          flags{f} = w;
          s = rmspace(s);
          [w s]=firstword(s);
          args{f} = w;
       end
     end
   end
  end
end 

%
% substitute the flags/args pairs in the .sc file
%
arg = firstarg;
while arg <= nargin
  eval(['Flag = var' int2str(arg) ';']);
  if arg == nargin
    fprintf('\nbinica(): Flag %s needs an argument.\n',Flag)
    return
  end
  eval(['Arg = var' int2str(arg+1) ';']);
  if strcmp(Flag,'pca') | strcmp(Flag,'PCA')
        ncomps = Arg; % get number of components out for reading wts.
  end
  arg = arg+2;

  nflags = f;
  for f=1:length(flags)   % replace arg with Arg
    if strcmp(Flag,flags{f})
       args{f} = num2str(Arg);
    end
  end
end

%
% insert name of data files, chans and frames
%
for x=1:length(flags)
  if strcmp(flags{x},'DataFile')
     args{x} = datafile;
  elseif strcmp(flags{x},'WeightsOutFile')
     weightsfile = ['binica' tmpint '.wts'];
     args{x} = weightsfile;
  elseif strcmp(flags{x},'SphereFile')
     spherefile = ['binica' tmpint '.sph'];
     args{x} = spherefile;
  elseif strcmp(flags{x},'chans')
     args{x} = int2str(nchans);
  elseif strcmp(flags{x},'frames')
     args{x} = int2str(nframes);
  end
end

%
% write the new .sc file
%
fid = fopen(scriptfile,'w');
for x=1:length(flags)
  fprintf(fid,'%s %s\n',flags{x},args{x});
end
fclose(fid);
if ~exist(scriptfile)
  fprintf('\nbinica(): ica script file %s not written.\n',...
                                   scriptfile);
  return
end
  
%
% %%%%%%%%%%%%%%%%%%%%%% run ica %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   fprintf('\nRunning ica from script file %s\n',scriptfile);
   eval(['!' ICA '<' scriptfile]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
% read in wts and sph results.
%
if ~exist('ncomps')
  ncomps = nchans;
end

wts = floatread(weightsfile,[ncomps Inf]);
if isempty(wts)
   fprintf('\nbinica(): weight matrix not read.\n')
   return
end
sph = floatread(spherefile,[nchans Inf]);
if isempty(sph)
   fprintf('\nbinica():  sphere matrix not read.\n')
   return
end
fprintf('\nica files left in pwd:\n');
eval(['!ls -l ' scriptfile ' ' weightsfile ' ' spherefile]);
fprintf('\n');

if isstr(data)
  whos wts sph
else
  whos data wts sph
end

%
% if created by binica(), rm temporary data file
%
if ~isempty(tmpdata)
  eval(['!rm ' datafile]);
end

%
%%%%%%%%%%%%%%%%%%% included functions %%%%%%%%%%%%%%%%%%%%%%
%
function sout = rmcomment(s,symb)
     n =1;
     while s(n)~=symb % discard # comments
        n = n+1;
     end
     if n == 1
        sout = [];
     else
       sout = s(1:n-1);
     end
    
function sout = rmspace(s)
       n=1;          % discard leading whitespace
       while n<length(s) & isspace(s(n))
          n = n+1;
       end
       if n<length(s)
          sout = s(n:end);
       else
          sout = [];
       end

function [word,sout] = firstword(s)
       n=1;         
       while n<=length(s) & ~isspace(s(n))
          n = n+1;
       end
       if n>length(s)
         word = [];
         sout = s;
       else
         word = s(1:n-1);
         sout = s(n:end);
       end
