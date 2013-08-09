% %this code does the MVDR beamforming for the data from GOM
% NoHydrophone = 64; 
% d = 1.5; f_samp = 8000; 
% sn = linspace(-1,1,401); %electrical angle (= sin theta) 
% ii = 10; 
% signal_ft = ones(18e3,64); 
% 
% delay = ([1:1:NoHydrophone]-32.5)*d*sn(ii); 
% delay_matrix = exp(-j*2*pi*[linspace(0, f_samp, size(signal_ft,1))]'*delay/1500);


function final_sum = MVDR(signal_ft, window1, aperture)


if strcmp(aperture, 'lf') || strcmp(aperture, 'LF')
    d = 1.5
elseif strcmp(aperture, 'hf') || strcmp(aperture, 'HF')
    d = 0.75; 
end

te = 0; 
sn = linspace(-1, 1, 401); 
N = length(sn);

f_samp = 8000; 
final_sum = zeros(N, size(signal_ft,1));
NoHydrophone = 64; 
T = size(signal_ft,1)/f_samp; 
freq = linspace(0,8000, size(rcv_t, 1)); 

window1 = lower(window1); 

if strcmp(window1, 'hanning') || strcmp(window1, 'hann') 
    window1 = hann(NoHydrophone); 
    window1 = window1/sum(window1); 
elseif strcmp(window, 'rect') || strcmp(window1, 'rectangular')
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
    s1 = delay; 
    delay_matrix = exp(-j*2*pi*[linspace(0, f_samp, size(signal_ft,1))]'*delay/1500);
    for jj = 1:size(signal_ft, 1)
        f = freq(jj); 
        s = exp(-j*2*pi*f*s1/1500); 
        s = transpose(s); 
        
        
        
    end
    
    
    rcv_spectrum2 = signal_ft.*delay_matrix;
    rcv_t_delay1 = fft(rcv_spectrum2)/T;     
    final_sum(te, :) = rcv_t_delay1*window1; 
    theta = asind(sn(ii)); 
end

%% estimating autocorrelation matrix 

for t_start = 0.1:1:(size(signal_ft,1)/8000 - 1)
    
end 

