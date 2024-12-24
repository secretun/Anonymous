%floatread() - Read matrix from float file.
%
%   >> A = floatread(filename,size,[format],[offset]) 
%
%   Read matrix a from specified
%   file while assuming four byte floating point numbers. The vector SIZE 
%   determine the number of float elements to be read and the dimensions
%   of the resulting matrix. If the last element of SIZE is INF the
%   size of the last dimension is determined by the file length. If
%   size is 'square,' floatread() attempts to read a square matrix.
%
%   The option FORMAT argument specifies the storage format as
%   defined by fopen(). Default format ([]) is 'native'.
%
%   The option OFFSET is offset in floats from the beginning of file (=0)
%   to start reading (4-byte floats assumed).
%   It uses fseek to read a portion of a large data file.
%
%   See also floatwrite(), fopen.

% xx-07-98  written by Sigurd Enghoff, Salk Institute
% 04-26-99  modified by Sigurd Enghoff to handle variable-sized and
%           multi-dimensional data.
% 07-08-99  modified by Sigurd Enghoff, FORMAT argument added.
% 02-08-00  help updated for toolbox inclusion -sm
% 02-14-00  added segment arg -sm
% 08-14-00  added size 'square' option -sm

function A = floatread(fname,Asize,fform,offset)

if nargin<2
  help floatread
  return
end

if ~exist('fform') | isempty(fform)|fform==0
	fform = 'native';
end

fid = fopen(fname,'rb',fform);
if fid>0 
 if exist('offset')
   stts = fseek(fid,4*offset,'bof');
   if stts ~= 0
     error('floatread fseek');
     return
   end
 end
 if ischar('Asize')
   if strcmp(Asize,'square')
         fseek(fid,0,'eof'); % go to end of file
         bytes = ftell(fid); % get byte position
         fseek(fid,0,'bof'); % rewind
         bytes = bytes/4; % nfloats
         froot = sqrt(bytes);
         if round(froot)*round(froot) ~= bytes
              error('floatread(): filelength is not square.')
         else
              Asize = [round(froot) round(froot)];
         end
   end
 end
 A = fread(fid,prod(Asize),'float');
else
 error('floatread() fopen');
 return
end
% fprintf('   %d floats read\n',prod(size(A)));

if Asize(end) == Inf
	Asize = Asize(1:end-1);
	A = reshape(A,[Asize length(A)/prod(Asize)]);
else
	A = reshape(A,Asize);
end

fclose(fid);
