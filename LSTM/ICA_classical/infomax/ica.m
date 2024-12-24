%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          MATLAB functions for psychophysiological data analysis           %%
%%                   including use of enhanced versions of the               %%
%%                 Independent Component Analysis (ICA) algorithm            %%
%%                         of Bell & Sejnowski (1995).                       %%
%%        By Scott Makeig, Colin Humphries, Sigurd Enghoff & Tzyy-Ping Jung, %%
%%        with contributions from Tony Bell, Sigurd Enghoff, Martin McKeown, %%
%%        Alex Dimitrov, Te-Won Lee, J-F Cardoso, Benjamin Blankertz et al.  %%
%%                        email: scott@salk.edu                              %%
%%               CNL/Salk Institute, 2000, Version 3.5.1                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  GENERAL ELECTROPHYSIOLOGICAL DATA PROCESSING TOOLS:
%
%  Make plot axes pop up into zoomable windows on mouse click      axcopy()
%  Find abs peak frames and amplitudes:                            abspeak()
%  Change reference from common to average:                        averef()
%  Simple block average data epochs:                               blockave()
%  Plot custom colorbar                                            cbar()
%  Make a 2-D scalp field movie:                                   eegmovie()
%  Frequency band filter data:                                     eegfilt()
%  View continuous data traces:                         eegplotg(),eegplot()
%  Average data epochs (with windowing options):                   erpave()
%  Display raw or smoothed single data epochs:          erpimage(),erpimages()
%  Re-align event-related epochs to given events:                  eventlock()
%  Plot one or more field maps on 3-D head model(s):    headplot(),compheads()
%  Construct a movie of a moving field on a 3-D head model:        headmovie()
%  Compute and view log power spectra of single data epochs:       logspec()
%  Select chans,frames,epochs of concatenated data epochs:         matsel()
%  Perform moving averaging on data:                               movav()
%  Plot a multichannel data epoch on a single axis:                ploterp()
%  Perform principal component analysis (PCA) via SVD              pcasvd()
%  Perform nonlinear (post-PCA) rotations:    varimax(), promax(), qrtimax()
%  View concatenated multichannel data epochs:                     plotdata()
%  View concatenated data epochs in topographic arrangement:       plottopo()
%  Plot a data epoch with topoplots at selected time points:       timtopo()
%  Change the data sampling rate:                                  resample()
%  Remove baseline means from data epochs:                         rmbase()
%  Regress out EEG data artifacts:                                 rmart()
%  View a 2-D or 3-D scalp-field movie:                            seemovie()
%  Time/frequency (ERSP, ITC) averages of single-trial data:       timef()
%  Iter-channel coherence averages of single-trial data:           crossf()
%  View data scalp topography(s):                       topoplot(),compmap()
%  View images using scalp topography info:                        imagetopo()
%  Convert Cartesian (x,y,z) channel locs to topoplot() format:    cart2topo()
%  Convert 2-D topoplot() channel locs to 3-D headplot() format:   topo2sph()
%  Convert 2-D headplot() channel locs to 2-D topoplot() format:   sph2topo()
%
%  SPECIFIC ICA TOOLS: 
%
%  Perform ICA analysis using logistic infomax or extended-infomax runica() 
%  Fast, compact Matlab MEX-file implementation of runica()        mexica() 
%  Fastest, most compact: system-call of binary runica()           binica()
%  Perform ICA analysis using 2nd & 4th-order cumulants (Cardoso)  jader()
%  Test ICA algorithm accuracy, varying data parameters:           testica()
%  Plot data and component envelopes:                    envproj() envtopo() 
%  Compute component activations:                                  icaact()
%  Compute component variances on scalp:                           icavar()
%  Make activations all rms-positive:                              posact()
%  Compute component projections:                                  icaproj() 
%  Plot the data decomposition:                      plotproj() -> chanproj()
%  Plot the data decomposition using plotopo():                    projtopo()
%  Sort ICA components by max projected latency and variance:      compsort()
%  Sort ICA components by mean projected variance only:            varsort() 
%  Compare ICA weight matrices:                       matcorr() -> matperm() 
%  Plot selected time periods of component activations:            tree()
%  View a projected ICA component (time course plus topo map):     compplot()
%  Squash or expand data into a PCA-defined subspace: pcsquash()-> expproj()
%                                                               -> pcexpand()
%  PACKAGE DEMO                                                    icademo()
%
%  REFERENCES:          http://www.cnl.salk.edu/~scott/icabib.html 
%  Further information: http://www.cnl.salk.edu/~scott/icafaq.html 
%  Please send news/bugs/fixes/suggestions to:  scott@salk.edu

help ica

% Begun: May 7, 1996
% Current version: Mon Aug 21 13:35:49 PDT 2000
