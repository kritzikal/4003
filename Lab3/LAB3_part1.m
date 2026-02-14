%% Pick on configuration and use these settings:
configMode = 2;  % 1=loopback, 2=two-dev/single-host, 3=TX host, 4=RX host

switch configMode
    case 1
        loopback = true;  runTx = true;  runRx = true;
    case 2
        loopback = false; runTx = true;  runRx = true;
    case 3
        loopback = false; runTx = true;  runRx = false;
    case 4
        loopback = false; runTx = false; runRx = true;
end

%% Creating RX and TX objects
% clear previous instances
clear tx rx

% add path to the directory where the function is found
%addpath('C:\Users\ebene\OneDrive\Desktop\PLUTO4Eben\common')

% Parameters
sampleRate = 30.72e6;
nsampsFrame = 2^12;

% Run the creation function
[tx, rx] = plutoCreateTxRx(createTx = runTx, createRx = runRx, loopback = loopback, ...
    nsampsFrame = nsampsFrame, sampleRate = sampleRate);

%% Create a Continuous Wave (CW) signal the length of one frame with a digital frequency
% of nu0 = 4/nsampFrame
nu0 = 4/nsampsFrame;

x = exp(2*pi*1i*nu0*(0:nsampsFrame-1)');


% Extract the real and imaginary parts
real_part = real(x);
imag_part = imag(x);

%% TODO - Plot real/imaginary of TX signal vs time
t_us = (0:nsampsFrame-1).' / sampleRate * 1e6;   % time in microseconds

figure;
subplot(2,1,1);
plot(t_us, real(x), 'LineWidth', 1);
grid on; xlabel('Time (\mus)'); ylabel('Amplitude');
title('TX Signal: Real Part');

subplot(2,1,2);
plot(t_us, imag(x), 'LineWidth', 1);
grid on; xlabel('Time (\mus)'); ylabel('Amplitude');
title('TX Signal: Imag Part');

%% Transmit continuously
if runTx
    tx.release();
    tx.transmitRepeat(x);
end

if ~runRx
    return;
end

%% Capturing Data Once
% We will now capture one frame of samples.  
% To capture the samples, we use the rx.capture() method.  
% This provides the samples in integer values.  
% We scale them to floating point and plot them.    
% If you are using two device mode, you should see a noisy version of the TX signal.

nbits = 12;               % number of ADC bits
fullScale = 2^(nbits-1);  % full scale
r = rx.capture(nsampsFrame);
r = single(r)/fullScale;

%% TODO - Plot real/imag of RX samples vs time
figure;
subplot(2,1,1);
plot(t_us, real(r), 'LineWidth', 1);
grid on; xlabel('Time (\mus)'); ylabel('Amplitude');
title('RX Samples: Real Part');

subplot(2,1,2);
plot(t_us, imag(r), 'LineWidth', 1);
grid on; xlabel('Time (\mus)'); ylabel('Amplitude');
title('RX Samples: Imag Part');


%% TODO - Print frame numbers where overflow occurred
nframes = 5;
data = complex(zeros(nsampsFrame, nframes, 'single'));
over = zeros(nframes,1);
for i = 1:nframes

    [r, valid, over(i)] = rx();
    if ~valid
        warning('Data is not valid');
    else
        data(:,i) = single(r)/ fullScale;
    end
    
end

idxOver = find(over ~= 0);

if isempty(idxOver)
    fprintf('No overflows detected.\n');
else
    fprintf('Overflow detected in frame(s): ');
    fprintf('%d ', idxOver);
    fprintf('\n');
end


%% Create dataConcat (length = nsampsFrame*nframes) & Plot real part(dataConcat) vs time(ms)
% Your data matrix is nsampsFrame x nframes. Concatenate column-wise:

dataConcat = data(:);

t_ms = (0:numel(dataConcat)-1).' / sampleRate * 1e3;   % time in milliseconds. Milliseconds help visualize longer-term continuity across frames

figure;
plot(t_ms, real(dataConcat), 'LineWidth', 1);
grid on;
xlabel('Time (ms)'); ylabel('Amplitude');
title('Real(dataConcat) vs Time (ms)');

%% Release the Tx/Rx Objects

release(tx); 
release(rx);
