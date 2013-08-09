%FUNCTION final_sum = steer_sperm_whale(signal_ft, window1, fs, noBeams, noHydrophones, sensor_positions, c, steer_angle)

%returns the beamtime data vs. angle 
%INPUT      -signal_ft: -a MxN complex Fourier transform of the received field on
%                        N phones, M is the total number of samples
%           -window   : -'rect' or 'hanning'. Default is Hanning window
%           -aperture : -'LF' for 415 and 'HF' for 950 aperture. Case-insensitive
%		-steer_angle: sin theta 
%OUTPUT     -final_sum: - MxnoBeams beamformed results corresponding to noBeams bearing angles equally spaced in sin(theta) from -1 to 1
%                           


function final_sum = steer_sperm_whale(signal_ft, window1, fs, noBeams, noHydrophones, sensor_positions, c, steer_angle)


% d = spacing; %hydrophone spacing 

% arrayCenter = mean(sensor_positions); 
arrayCenter = (sensor_positions(end) + sensor_positions(1))/2 + sensor_positions(1); 

 
final_sum = zeros(1, size(signal_ft,1)); %initialize an array of steered beams
T = size(signal_ft,1)/fs; %total signal duration 

window1 = lower(window1); 

if strcmp(window1, 'hanning') || strcmp(window1, 'hann') 
    window1 = hann(noHydrophones); 
    window1 = window1/sum(window1); 
elseif strcmp(window1, 'rect') || strcmp(window1, 'rectangular')
    window1 = ones(noHydrophones, 1); 
    window1 = window1/sum(window1); 
end
display(['steering to sin theta = ' num2str(steer_angle)])
% total_energy = []; 
% max1 = -10000000; 

delay = (sensor_positions - arrayCenter)*steer_angle;
delay_matrix = exp(-j*2*pi*[linspace(0, fs, size(signal_ft,1))]'*delay/c);
rcv_ft2 = signal_ft.*delay_matrix; %delayed signals Fourier transform matrix
rcv_t_delay1 = fft(rcv_ft2)/T;   %delayed signals in time  
final_sum = rcv_t_delay1*window1; %sum after delaying 

display('done'); 
