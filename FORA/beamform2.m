%FUNCTION final_sum = beamform2(signal_ft, window, aperture)
%returns the beamtime data vs. angle 
%INPUT      -signal_ft: -a Mx64 complex Fourier transform of the received field on
%                        64 phones, M is the total number of samples
%           -window   : -'rect' or 'hanning'. Default is Hanning window
%           -aperture : -'LF' for 415 and 'HF' for 950 aperture. Case-insensitive
%
%OUTPUT     -final_sum: - Mx401, 401 beamformed results corresponding to 401 bearing angles equally spaced in sin(theta) from -1 to 1
%                           
% Modified Jan 24 '11: eliminate inefficient loops and speed up by about 1.7 times.
% Thanks to Dave for pointing this out. 
% last modified: Jun 28 '2011 by DD. Minor fixes. 

function final_sum = beamform2(signal_ft, window1, aperture)


if strcmp(aperture, 'lf') || strcmp(aperture, 'LF')
    d = 1.5;
elseif strcmp(aperture, 'mf') || strcmp(aperture, 'MF')
    d = 0.75; 
elseif strcmp(aperture, 'hf') || strcmp(aperture, 'HF')
    d = 0.375; 
end

te = 0; 
sn = linspace(-1, 1, 401); 
N = length(sn);

f_samp = 8000; 
final_sum = zeros(N, size(signal_ft,1));
NoHydrophone = 64; 
T = size(signal_ft,1)/f_samp; 

window1 = lower(window1); 

if strcmp(window1, 'hanning') || strcmp(window1, 'hann') 
    window1 = hann(NoHydrophone); 
    window1 = window1/sum(window1); 
elseif strcmp(window1, 'rect') || strcmp(window1, 'rectangular')
    window1 = ones(NoHydrophone, 1); 
    window1 = window1/sum(window1); 
end

total_energy = []; 
max1 = -10000000; 
for ii = 1:length(sn)
    te = te +1;    
    display('beamforming')
    display([num2str(ii/length(sn)*100) , '% complete']);   
    delay = ([1:1:NoHydrophone]-32.5)*d*sn(ii); 
    delay_matrix = exp(-j*2*pi*[linspace(0, f_samp, size(signal_ft,1))]'*delay/1500);
    rcv_spectrum2 = signal_ft.*delay_matrix;
    rcv_t_delay1 = fft(rcv_spectrum2)/T;     
    final_sum(te, :) = rcv_t_delay1*window1; 
    theta = asind(sn(ii)); 
end