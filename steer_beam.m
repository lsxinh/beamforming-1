function final_sum = steer_beam(signal_ft, filter_window, aperture, steer_angle)
%FUNCTION final_sum = steer_beam(signal_ft, filter_window, aperture, steer_angle)
% is a modified version of beamform2 when the beam angle has been
% accurately known, and the beamformed signal over a long duration is
% needed. This is to save the time of beamforming over all 360 degrees. 
%INPUT      -signal_ft: -a Mx64 complex Fourier transform of the received field on
%                        64 phones, M is the total number of samples
%           -window   : -'rect' or 'hanning'. Default is Hanning window
%           -aperture : -'LF' for 415 and 'MF' for 950 aperture, HF for high frequency aperture
%                       . Case-insensitive
%           -steer_angle: sine value of the actual angle. eg. 0.5 (instead of 30 degrees)
%
%OUTPUT     -final_sum: - Mx1,  beamformed results corresponding to the
%bearing angle specified 

% Modified Jan 24 '11: eliminate inefficient loops and speed up by about 1.7 times.
% Thanks to Dave for pointing this out. 
% last modified: Jun28 '2012 by DD. Minor fixes.


if strcmp(aperture, 'lf') || strcmp(aperture, 'LF')
    d = 1.5
elseif strcmp(aperture, 'mf') || strcmp(aperture, 'mF')
    d = 0.75; 
elseif strcmp(aperture, 'hf') || strcmp(aperture, 'HF')
    d = 0.375; 
end

te = 0; 
sn = steer_angle; 

f_samp = 8000; 
final_sum = zeros(1, size(signal_ft,1));
NoHydrophone = 64; 
T = size(signal_ft,1)/f_samp; 

if strcmp(filter_window, 'Hanning') || strcmp(filter_window, 'hanning') || strcmp(filter_window, 'hann')
    window1 = hann(NoHydrophone); 
    window1 = window1/sum(window1); 
elseif strcmp(filter_window, 'rect') || strcmp(filter_window, 'rectangular')
    window1 = ones(NoHydrophone, 1); 
    window1 = window1/sum(window1); 
end

% display('beamforming')
% display([num2str(ii/length(sn)*100) , '% complete']);   
delay = ([1:1:NoHydrophone]-32.5)*d*sn; 
delay_matrix = exp(-j*2*pi*[linspace(0, f_samp, size(signal_ft,1))]'*delay/1500);
rcv_spectrum2 = signal_ft.*delay_matrix; %multiplication with exp(j2pif dt) in freq domain corresponds to shift in time 
rcv_t_delay1 = fft(rcv_spectrum2)/T;     
final_sum = rcv_t_delay1*window1; 



