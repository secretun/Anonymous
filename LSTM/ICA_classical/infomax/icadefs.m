%
% icadefs - file of default filenames and constants to source in the ICA /ERP
%               package functions.  Insert local dir reference below. 

% 05-20-97   Scott Makeig, CNL / Salk Institute, La Jolla CA
% 08-04-00  added ICA and SC -sm

if isunix
 ICADIR = [ '/home/scott/matlab/' ];   % INSERT Unix Matlab ICA dirname here
else
 ICADIR = [ 'd:\important\ica51\2001\ica5.3.6\' ];      % INSERT PC matlab ICA dirname here
end

ICA = '/home/enghoff/stalone/sica';    % INSERT ica executable for binica.m

SC  =  [ICADIR 'binica.sc'];           % master .sc file for binica.m

ENVCOLORS  = [ ICADIR 'envproj.col' ]; % default color-order
%                                        filename for envproj.m here.

PROJCOLORS = [ ICADIR 'white1st.col' ];% default color-order
%                                         filename for plotproj.m here.
BACKCOLOR  = [0.7 0.7 0.7];            % background color for plotting

MAXENVPLOTCHANS   = 256;  % maximum number of channels to plot in envproj.m
MAXPLOTDATACHANS  = 256;  % maximum number of channels to plot in dataplot.m
                          %         and functions that use it.
MAXPLOTDATAEPOCHS = 256;  % maximum number of epochs to plot in dataplot.m
MAXEEGPLOTCHANS   = 256;  % maximum number of channels to plot in eegplot.m
MAXTOPOPLOTCHANS  = 256;  % maximum number of channels to plot in topoplot.m
DEFAULT_ELOC  = 'chan_file'; % default electrode location file for topoplot.m
DEFAULT_SRATE = 256;      % default sampling rate
DEFAULT_EPOCH = 10;       % default max epoch width to plot in eegplot(s).m


if strcmp(ICADIR,'XXX/')
    fprintf('===============================================\n');
    fprintf('You have not set the ICADIR variable in icadefs.m\n');
    fprintf('===============================================\n');
    return
end

