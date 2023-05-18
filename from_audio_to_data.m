function [x,yclean,fs,N] = from_audio_to_data(File,y_max_length,fs_)

% Downsample the signal in 'File' to the frenquency 'fs_' 

[y,fs] = audioread([File,'.wav']); % reads in the file
yclean = resample(y, fs_, fs); % downsample the input
fs = fs_; % reducing the sampling frequency
N=length(yclean); % signal size
x=linspace(0,y_max_length,N); % time instances


end

