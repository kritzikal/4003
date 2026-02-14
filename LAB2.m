%% Lab 2 — Fourier Analysis, FFT, Sampling Theory, and Convolution
% Course: Wireless Communications
% Duration: ~2 hours 
% Coverage mapping (textbooks):
% • Principles of Modern Communication Systems — Sec. 2.5–2.8 (Fourier analysis), Sec. 2.8.9 (Convolution)
% • SDR for Engineers — Sec. 2.1.3 (FFT), Sec. 2.2 & 5.2 (Sampling Theory)

%% 0) Setup (Run Once)
clear; close all; clc;
fs  = 50e3;                 % base sample rate (50 kHz) for most demos
Ts  = 1/fs;
T   = 0.2;                  % short duration to keep figures crisp
t   = 0:Ts:T-Ts;
rng(7);

%% 1) Fourier Analysis — Fourier Series & Reconstruction (Principles 2.5–2.7)
% Goal: Reconstruct periodic signals from sinusoids; inspect convergence and Gibbs effect;
% verify Parseval’s relation; observe faster convergence for smoother waveforms.

% 1.1 Square Wave: Synthesis by Odd Harmonics (Gibbs)
f0 = 50;                % fundamental frequency (Hz)
w0 = 2*pi*f0;
t1 = 0:Ts:0.1-Ts;       % 100 ms window
x_true = square(w0*t1); % ideal ±1 square wave

% Partial sums: N harmonics (odd only)
Ns = [1, 3, 5, 9, 25, 101];
figure('Name','Square Wave Fourier Series (Partial Sums)');
for i = 1:numel(Ns)
    Nhar = Ns(i);
    x_fs = zeros(size(t1));
    for k = 1:Nhar
        n = 2*k-1; % odd
        x_fs = x_fs + (4/pi)*(1/n)*sin(n*w0*t1);
    end
    subplot(numel(Ns),1,i);
    h1 = plot(t1, x_true, ':', 'LineWidth', 1.2); hold on;
    h2 = plot(t1, x_fs,  '--', 'LineWidth', 1.2); grid on;
    xlim([t1(1) t1(end)]);
    ylim([-1.5, 1.5]);
    title(sprintf('Odd Harmonics: N=%d (up to %d f0)', Nhar, 2*Nhar-1));
    
    if i == 1
        legend([h1 h2], {'True square','Fourier partial sum'}, ...
               'Location','northeast');
    end
    
    if i==numel(Ns), xlabel('Time (s)'); end
end


%% Next -- 1.2 Triangle Wave:
x_tri = sawtooth(w0*t1, 0.5);     % triangle in [-1,1]

Ns2 = [1, 3, 5, 9, 25];
figure('Name','Triangle Wave Fourier Series');
for i = 1:numel(Ns2)
    Nhar = Ns2(i);
    x_fs = zeros(size(t1));
    for k = 1:Nhar
        n = 2*k-1; % odd
        x_fs = x_fs + (8/pi^2)*(1/n^2)*(-1)^((n-1)/2)*sin(n*w0*t1);
    end
    subplot(numel(Ns2),1,i);
    h1 = plot(t1, x_tri, ':', 'LineWidth', 1.2); hold on;
    h2 = plot(t1, x_fs,  '--', 'LineWidth', 1.2); grid on;
    xlim([t1(1) t1(end)]);
    title(sprintf('Triangle: N=%d odd harmonics (1/n^2 decay)', Nhar));
    ylim([-1.5, 1.5]);
    
    if i == 1
        legend([h1 h2], {'True triangle','Fourier partial sum'}, ...
               'Location','northeast');
    end
    
    if i==numel(Ns2), xlabel('Time (s)'); end

end

%% Next -- 1.3 Parseval Check (average power over window)
Nhar = 25;
x_fs = zeros(size(t1));
for k = 1:Nhar
    n = 2*k-1;
    x_fs = x_fs + (4/pi)*(1/n)*sin(n*w0*t1);
end

Tobs   = numel(t1)*Ts;
P_time = (1/Tobs) * trapz(t1, x_fs.^2);

amps   = (4/pi)./(1:2:(2*Nhar-1));
P_coeff = 0.5*sum(amps.^2);

fprintf('Parseval check (square, N=%d): Time-domain P ≈ %.4f, Coeff-domain P ≈ %.4f\n', ...
    Nhar, P_time, P_coeff);

%% Part 2) Fourier Transform & Spectral Analysis (Principles 2.7; SDR4E 2.1.3)
% Compare spectra for different window lengths and windows; interpret bin spacing,
% leakage, and zero-padding.

%% 2.1a Two-Tone Test & Window Effects
f1 = 700; f2 = 1200; A1 = 1; A2 = 0.6; phi2 = pi/4;
x = A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t + phi2);

%Case A: window = rect (no window), Nfft = length(x)
NfftA = numel(x);
[fA, AA, ~] = ampSpec(x, fs);

%Case B: apply Hann window (reduce leakage)
w = hann(numel(x))'; xw = x.*w;
[fB, AB, ~] = ampSpec(xw, fs);

% Plot single-sided spectra up to 5 kHz
limF = 5000;
figure('Name','Windowing & Leakage');
subplot(2,1,1);
plot(fA, AA); xlim([0, limF]); grid on; title('Rectangular Window');
ylabel('|X(f)| (linear)');
subplot(2,1,2);
plot(fB, AB); xlim([0, limF]); grid on; title('Hann Window');
xlabel('Frequency (Hz)'); ylabel('|X(f)| (linear)');


%% 2.1b Two-Tone Test & Window Effects
f1 = 733; f2 = 1217; A1 = 1; A2 = 0.6; phi2 = pi/4;   % intentionally non-bin-aligned
x  = A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t + phi2);

% Rectangular (no window)
[fA, AA, ~] = ampSpec(x, fs);

% Hann window (reduced leakage) + simple amplitude compensation
w  = hann(numel(x))';
xw = (x .* w) / mean(w);
[fB, AB, ~] = ampSpec(xw, fs);

% Plot single-sided spectra up to 5 kHz in dB
limF = 5000;
figure('Name','Windowing & Leakage');
subplot(2,1,1);
plot(fA, 20*log10(AA + eps)); xlim([0 limF]); grid on;
title('Rectangular Window'); ylabel('|X(f)| (dB)'); ylim([-120 10]);

subplot(2,1,2);
plot(fB, 20*log10(AB + eps)); xlim([0 limF]); grid on;
title('Hann Window'); xlabel('Frequency (Hz)'); ylabel('|X(f)| (dB)'); ylim([-120 10]);


%% 2.2 Phase Retrieval & Amplitude Calibration
% Recover amplitude/phase of a single tone with perfect bin alignment
f0 = 1000; Nper = round(fs/f0); Nrec = 4*Nper; t2 = (0:Nrec-1)/fs;
x2 = 1.5*cos(2*pi*f0*t2 + 0.3);
X2 = fft(x2); k0 = round(f0*Nrec/fs)+1;   % +1 for MATLAB indexing
A_est = 2*abs(X2(k0))/Nrec;               % single-sided amplitude
phi_est = angle(X2(k0));
fprintf('Tone at %d Hz: true A=1.5, est A=%.3f; true φ=0.300, est φ=%.3f\n', f0, A_est, phi_est);

%% 3) Fast Fourier Transform (SDR4Engineers 2.1.3)
% Practical FFT use: fft/ifft, fftshift, circular vs linear convolution, overlap-add outline.

%% 3.1 fft, ifft, and fftshift Basics
x3 = cos(2*pi*800*t) + 0.5*sin(2*pi*1400*t);
X3 = fft(x3);
X3s= fftshift(X3);
[f2, ~] = freqAxis(numel(x3), fs);
figure('Name','fft/fftshift Basics');
plot(f2, abs(X3s)/numel(x3)); grid on; xlim([-4000, 4000]);
xlabel('Frequency (Hz)'); ylabel('Amplitude (scaled)'); title('Two-Sided |X(f)| (fftshifted)');

%% 4) Sampling Theory (SDR4Engineers 2.2, 5.2)
% Nyquist condition, aliasing, anti-aliasing filters, interpolation/reconstruction intuition.

%% 4.1 Aliasing Demonstration (fs high vs near Nyquist vs low)
f_sig = 9e3;                   % 9 kHz tone
Tdemo = 0.01;                  % 10 ms demo
ts = 0:1/200e3:Tdemo;          % "continuous" reference at 200 kHz
x_ct = cos(2*pi*f_sig*ts);

% High fs (no alias), near Nyquist, and below Nyquist
fsA = 50e3;  nA = 0:round(fsA*Tdemo)-1; tA = nA/fsA; xA = cos(2*pi*f_sig*tA);
fsB = 18e3;  nB = 0:round(fsB*Tdemo)-1; tB = nB/fsB; xB = cos(2*pi*f_sig*tB);
fsC = 12e3;  nC = 0:round(fsC*Tdemo)-1; tC = nC/fsC; xC = cos(2*pi*f_sig*tC);

figure('Name','Aliasing Demo');
subplot(3,1,1); plot(ts, x_ct, 'LineWidth', 1); hold on; stem(tA, xA, '.'); hold off;
title('fs = 50 kHz (well above Nyquist)'); grid on;
subplot(3,1,2); plot(ts, x_ct, 'LineWidth', 1); hold on; stem(tB, xB, '.'); hold off;
title('fs = 18 kHz (near Nyquist 2*9 kHz)'); grid on;
subplot(3,1,3); plot(ts, x_ct, 'LineWidth', 1); hold on; stem(tC, xC, '.'); hold off;
title('fs = 12 kHz (below Nyquist) → alias'); grid on; xlabel('Time (s)');

%% Next -- 4.2 Anti-Aliasing Filter before Downsampling
fs_high = 100e3; tH = 0:1/fs_high:Tdemo-1/fs_high;
xH = cos(2*pi*9e3*tH) + 0.5*cos(2*pi*22e3*tH);   % include an out-of-band tone
fs_low = 24e3; R = round(fs_high/fs_low);        % integer factor ~ 4.17 ⇒ use exact 4 for demo
R = 4; fs_low = fs_high/R;
% Design LPF with cutoff < fs_low/2 to satisfy Nyquist at low rate
hAA = fir1(128, (fs_low/2 - 1000)/(fs_high/2));  % margin of 1 kHz
xF  = filtfilt(hAA, 1, xH);
xD  = xF(1:R:end); tD = tH(1:R:end);

figure('Name','Anti-Aliasing before Decimation');
subplot(2,1,1);
pwelch(xH, [], [], [], fs_high); title('Before AA filter');
subplot(2,1,2);
pwelch(xD, [], [], [], fs_low); title('After AA + Decimate');

%% Next -- 4.3 Interpolation / Reconstruction (sinc intuition)
% Reconstruct a lowpass signal sampled at fsA using a zero-order hold vs. ideal (sinc) LPF.
fsA = 40e3; f_lp = 4e3; Trec = 0.02;
tA = 0:1/fsA:Trec-1/fsA; x_lp = lowpass(cos(2*pi*3e3*tA)+0.7*cos(2*pi*1e3*tA+0.5), f_lp, fsA);

% Zero-order hold (staircase) for intuition
tHR = 0:1/200e3:Trec-1/200e3;
x_zoh = interp1(tA, x_lp, tHR, 'previous'); % ZOH visualization
% Ideal-like LPF reconstruction via upsampling & FIR lowpass
L = 5; x_up = upsample(x_lp, L);
hL = fir1(128, (f_lp)/(L*fsA/2));  % cutoff at f_lp
x_rec = filtfilt(hL, 1, x_up);     % simple two-sided filtering
tREC = (0:numel(x_rec)-1)/(L*fsA);

figure('Name','Reconstruction Intuition');
subplot(2,1,1); plot(tHR, x_zoh); grid on; title('ZOH (qualitative)');
subplot(2,1,2); plot(tREC, x_rec); grid on; title('LPF Reconstructed (upsample+filter)'); xlabel('Time (s)');

%% 5) Convolution (Principles 2.8.9)
% Time-domain convolution, convolution theorem, linear vs circular, and a simple LTI example.

%% 5.1 RC-like Impulse Response and Convolution
RC = 1e-3;                         % 1 ms
t5 = 0:Ts:0.02-Ts;
h_rc = (1/RC)*exp(-t5/RC);         % continuous-time analog shape sampled
x5 = square(2*pi*500*t5);          % input square at 500 Hz
y5 = conv(x5, h_rc)*Ts;            % approximate continuous-time conv (∫ x(τ)h(t-τ)dτ)
t5y = (0:numel(y5)-1)*Ts;

figure('Name','RC Response to Square');
plot(t5, x5, '--'); hold on; plot(t5y, y5); grid on;
legend('Input (square)','Output (conv with h_{RC})'); xlabel('Time (s)'); title('Low-pass smoothing by RC');


%% Helper utilities (used throughout)
% freqAxis(): build two-sided and one-sided frequency axes
% ampSpec(): amplitude-correct single-sided amplitude spectrum for real signals
% zpad(): zero-pad to next power of 2 (or specified length)
function [f2, f1] = freqAxis(N, fs)
    % Two-sided (fftshifted) and one-sided axes in Hz
    f2 = (-N/2:N/2-1)*(fs/N);
    f1 = (0:N/2)*(fs/N);
end

function [f1, A1, P1] = ampSpec(x, fs)
    % Single-sided amplitude spectrum with amplitude correction for real x
    N  = numel(x);
    X  = fft(x);

    A2 = abs(X)/N;                  % two-sided amplitude
    P2 = angle(X);
    A1 = A2(1:N/2+1);
    P1 = P2(1:N/2+1);
    A1(2:end-1) = 2*A1(2:end-1);    % fold negative freqs
    f1 = (0:N/2)*(fs/N);
end

function y = zpad(x, N)
    if nargin<2, N = 2^nextpow2(numel(x)); end
    y = zeros(1,N); y(1:numel(x)) = x;
end

%% 6) Suggested Analytical Questions (by Topic)
% Fourier Series / Fourier Transform
%  FQ1: For a periodic signal with Fourier coefficients a_n, b_n, derive the average power from coefficients.
%       Answer: P = (a_0^2)/4 + (1/2)∑_{n=1}^∞ (a_n^2 + b_n^2).
%  FQ2: Show that for a square wave, odd harmonics dominate. Why are even harmonics zero?
%       Answer: Half-wave symmetry ⇒ even harmonics vanish by orthogonality.
%  FQ3: Given a Hann window, what is the mainlobe width (in bins) and sidelobe attenuation (approximate)?
%       Answer: ~2 bins mainlobe (full width between first nulls ≈ 4π/N rad), sidelobes ≈ −31 dB.

% FFT
%  FFTQ1: Prove complexity difference between DFT (O(N^2)) and FFT (O(N log N)).
%         Answer: Divide-and-conquer splits into sub-DFTs of size N/2 recursively; sum work across stages.
%  FFTQ2: Explain why bin centers occur at k*fs/N and how to align a tone exactly to a bin.
%         Answer: Make the record contain an integer number of periods so f0 = m*fs/N.

% Sampling Theory
%  SQ1: Derive aliasing mapping formula f_alias = |f − k*fs| restricted to [0, fs/2].
%       Answer: From periodic replication of spectra at multiples of fs and folding into baseband.
%  SQ2: Show that ideal reconstruction filter has frequency response H(f) = rect(f/fs).
%       Answer: Sinc in time ↔ rectangular LPF in frequency; passband up to fs/2.

% Convolution
%  CQ1: Show that linear convolution length is N+M−1 for sequences of lengths N and M.
%       Answer: Support of sums over shifts extends from 0 to (N−1)+(M−1).
%  CQ2: For an FIR with linear phase, show group delay (N−1)/2.
%       Answer: Symmetric coefficients ⇒ constant phase slope.

% End of Lab 2
