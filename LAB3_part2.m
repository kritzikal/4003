%% Lab 3 — Intro to SDR with ADALM‑PLUTO: Spectrum Analyzer
% Course: Wireless Communications
% Duration: ~2 hours

%% 0) Setup & Device Check (Run Once)
clear; close all; clc;
useHardware = true;

fs   = 1.0e6;
fc   = 900e6;
gain = 50;

devFound = false;
devURI   = '';

if useHardware
    try
        radios = findPlutoRadio();            % <-- struct array

        if isempty(radios)
            error('findPlutoRadio returned empty (no Pluto detected).');
        end

        devURI   = radios(1).RadioID;         % <-- e.g. 'usb:0'
        devFound = true;

        fprintf('Pluto found at %s\n', devURI);

    catch ME
        warning('Pluto not detected: %s\nFalling back to synthetic signal.', ME.message);
        useHardware = false;
    end
end
%% 1) Create Receiver & Spectrum Analyzer
% Build an sdrrx('Pluto') object (hardware) or a synthetic source (software).

if useHardware
    rx = sdrrx('Pluto', ...
        'RadioID',            devURI, ...
        'CenterFrequency',    fc, ...
        'BasebandSampleRate', fs, ...
        'SamplesPerFrame',    4096, ...
        'GainSource',         'Manual', ...
        'Gain',               gain);
else
    % Synthetic RF: mix a few passband tones around fc (simulated downconversion already applied)
    rx = [];  % placeholder
end

% Create Spectrum Analyzer
spec = spectrumAnalyzer( ...
    'Name','Pluto Spectrum', ...
    'SampleRate',fs, ...
    'AveragingMethod','exponential', ...
    'RBWSource','auto', ...
    'PlotAsTwoSidedSpectrum', true, ...
    'FrequencySpan','Full', ...
    'YLimits', [-30 -20], ...
    'SpectrumType','Power density');  % dBW/Hz by default


%% 2) Live Spectrum — Tune, Gain, & Sample Rate Experiments
% Run for a short window; adjust fc, fs, gain; observe spectrum.

runSeconds = 15;     % shorten for demo; increase in class
frames = ceil(runSeconds * fs / 4096);

fprintf('Running live spectrum... (%.0f s)\n', runSeconds);
for k = 1:frames
    if useHardware
        x = rx();                            % complex baseband samples
    else
        % Synthetic: two tones + noise in baseband to emulate nearby stations
        N = 4096;
        t = (0:N-1).'/fs;
        x = 0.5*exp(1j*2*pi*100e3*t) + 0.2*exp(-1j*2*pi*180e3*t) + 0.05*(randn(N,1)+1j*randn(N,1));
    end
    spec(x);
end
release(spec);

%% 3) DC Spur Mitigation (Center Spike)
% Many direct-conversion radios exhibit a DC spike at 0 Hz (center). Demonstrate simple mitigation.

fullScale = 2^(12-1);   % Pluto ADC signed 12-bit => 2048

if useHardware
    x = rx();                           % complex int16
    x = single(x) / fullScale;          % convert + scale to float
else
    N = 4096; t = (0:N-1).'/fs;
    x = 0.5*exp(1j*2*pi*10e3*t) + 0.02*(randn(N,1)+1j*randn(N,1));
end

x_dc = x - mean(x);                      % remove DC offset (works now)

spec2 = spectrumAnalyzer( ...
    'Name','DC Spur Mitigation', ...
    'SampleRate', fs, ...
    'SpectralAverages', 16, ...
    'YLimits', [-120 -20]);              % typical PSD view

spec2(x); 
spec2(x_dc);
release(spec2);

%% 4) Waterfall & Frequency Sweep (Band Survey)
fStart = 902e6; fStop = 928e6; nSteps = 41;
fc_list = linspace(fStart, fStop, nSteps);

Nwelch  = 4096; 
nover   = 2048;
w       = hann(Nwelch,'periodic');

% For COMPLEX IQ with 'centered', PSD length is Nwelch (NOT Nwelch/2+1)
waterfallPSD = zeros(nSteps, Nwelch);

% Frequency axis for centered PSD: [-fs/2, fs/2)
fAxis = (-Nwelch/2:Nwelch/2-1)*(fs/Nwelch);   % 1 x Nwelch

for i = 1:nSteps
    if useHardware
        rx.CenterFrequency = fc_list(i);
        pause(0.05);
        x = rx();                             % int16 complex IQ from Pluto
        x = single(x) / 2^(12-1);             % scale to ~[-1,1], Pluto ADC is 12-bit signed
    else
        N = 4096; t = (0:N-1).'/fs;
        f_offset = (98.1e6 - fc_list(i));     
        x = 0.5*exp(1j*2*pi*f_offset*t) + 0.02*(randn(N,1)+1j*randn(N,1));
    end

    % Centered PSD for complex signals
    [Pxx, ~] = pwelch(x, w, nover, Nwelch, fs, 'centered');   % Nwelch x 1
    waterfallPSD(i,:) = 10*log10(Pxx.' + eps);                % 1 x Nwelch
end

% Absolute RF axis: fc + baseband offset
absFreq = fc_list(:) + fAxis;   % nSteps x Nwelch (implicit expansion)

figure('Name','Waterfall Sweep (FM Band)');
imagesc(absFreq(1,:)/1e6, 1:nSteps, waterfallPSD); axis xy;
xlabel('RF Frequency (MHz)'); ylabel('Sweep Index'); title('Welch PSD Waterfall (Centered)');
colorbar;
%colormap turbo;

%% 5) Channel Power & Occupancy (Quick Measurements)
if useHardware
    x = rx();
    x = single(x) / 2^(12-1);
else
    N=4096; t=(0:N-1).'/fs;
    x = 0.6*exp(1j*2*pi*30e3*t) + 0.02*(randn(N,1)+1j*randn(N,1));
end

% Use centered PSD for complex IQ (recommended)
Nfft = 4096;
[Pxx, fvec] = pwelch(x, hann(Nfft,'periodic'), Nfft/2, Nfft, fs, 'centered'); % Pxx in W/Hz

% Find strongest peak
[~, kmax] = max(Pxx);

BW = 200e3;                         % channel width
f0 = fvec(kmax);                    % peak frequency (Hz) in baseband

% Indices covering [f0-BW/2, f0+BW/2]
fL = find(fvec >= (f0 - BW/2), 1, 'first');
fH = find(fvec <= (f0 + BW/2), 1, 'last');

% Channel power by integrating PSD over frequency
df = fvec(2) - fvec(1);
chanPower = sum(Pxx(fL:fH)) * df;    % watts
chanPower_dBW = 10*log10(chanPower + eps);

fprintf('Estimated channel power over ~%.0f kHz around peak: %.1f dBW\n', BW/1e3, chanPower_dBW);

%% 6) Record IQ to File (for Offline Analysis)
captureSeconds = 2; frames = ceil(captureSeconds * fs / 4096);
Xcap = complex([]);
for k = 1:frames
    if useHardware, x = rx(); else, N=4096; t=(0:N-1).'/fs; x = 0.4*exp(1j*2*pi*50e3*t)+0.02*(randn(N,1)+1j*randn(N,1)); end
    Xcap = [Xcap; x]; %#ok<AGROW>
end
save('pluto_capture.mat','Xcap','fs','fc','gain');
fprintf('Saved capture to pluto_capture.mat (%.1f s, %.0f samples)\n', numel(Xcap)/fs, numel(Xcap));

%% 7) Clean Up
if useHardware
    release(rx);
end
disp('Lab 3 complete.');


% — End of Lab 3 —
