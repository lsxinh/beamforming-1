%FUNCTION final_sum = beamform(signal_ft, window1, fs, noBeams, noHydrophones, sensor_positions, c)
%returns the beamtime data vs. angle 
%INPUT      -signal_ft: -a MxN complex Fourier transform of the received field on
%                        N phones, M is the total number of samples
%           -window1   : -'rect' or 'hanning'. Default is Hanning window
%	    -fs : sampling frequency
%           -noBeams: number of desired beams. At least = no. of hydrophones. 
%	    -noHydrophones = number of hydrophones
%	    -sensor_position: positions of array sensors. 
%           -aperture : -'LF' for 415 and 'HF' for 950 aperture. Case-insensitive
%
%OUTPUT     -final_sum: - M x noBeams, beamformed results corresponding to 401 bearing angles equally spaced in sin(theta) from -1 to 1
%                           
% Last updated by DD Tran, Aug 8, 2013.  

function final_sum = beamform(signal_ft, window1, fs, noBeams, noHydrophones, sensor_positions, c)


% d = spacing; %hydrophone spacing 

% arrayCenter = mean(sensor_positions); 
arrayCenter = (sensor_positions(end) + sensor_positions(1))/2 + sensor_positions(1); 

sn = linspace(-1, 1, noBeams); %steered directions
 
final_sum = zeros(noBeams, size(signal_ft,1)); %initialize an array of steered beams
T = size(signal_ft,1)/fs; %total signal duration 

window1 = lower(window1); 

if strcmp(window1, 'hanning') || strcmp(window1, 'hann') 
    window1 = hann(noHydrophones); 
    window1 = window1/sum(window1); 
elseif strcmp(window1, 'rect') || strcmp(window1, 'rectangular')
    window1 = ones(noHydrophones, 1); 
    window1 = window1/sum(window1); 
end
display('beamforming')
% total_energy = []; 
% max1 = -10000000; 
te = 0; %counter for beam index
for ii = 1:length(sn)
    te = te +1;    
   
%     display([num2str(ii/length(sn)*100) , '% complete']);   
%     delay = ([1:1:noHydrophones]-(noHydrophones + 1)/2)*d*sn(ii); 
    delay = (sensor_positions - arrayCenter)*sn(ii);
    delay_matrix = exp(-j*2*pi*[linspace(0, fs, size(signal_ft,1))]'*delay/c);
    rcv_ft2 = signal_ft.*delay_matrix; %delayed signals Fourier transform matrix
    rcv_t_delay1 = fft(rcv_ft2)/T;   %delayed signals in time  
    final_sum(te, :) = rcv_t_delay1*window1; %sum after delaying 
%     theta = asind(sn(ii)); 
end
display('done'); 
