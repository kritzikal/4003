%% Lab 1 — Communication Basics & Signals
% Course: Wireless Communications
% Duration: ~2 hours

%% 1) Pre-Lab Setup (Run Once)
clear; close all; clc;
fs  = 1e5;                 % base sample rate (100 kHz) for baseband demos
Ts  = 1/fs;
T   = 1;                   % default 1 s duration
t   = (0:fs*T-1)*Ts;
rng(42);                   % seed for reproducibility

%% 2) Classification of Signals (Ch. 2.1–2.3)
% 2.1 Continuous-Time vs Discrete-Time

f0 = 500;                          % 500 Hz tone (continuous-time)

tc = linspace(0, 0.02, 10e4);      % 'continuous-time' grid (very fine sampling)

xc = cos(2*pi*f0*tc);              % Continuous Time reference cosine waveform

fs1 = 8e3;  n1 = 0:round(fs1*0.02)-1; td1 = n1/fs1; xd1 = cos(2*pi*f0*td1);
fs2 = 600;  n2 = 0:round(fs2*0.02)-1; td2 = n2/fs2; xd2 = cos(2*pi*f0*td2);

figure('Name','CT vs DT'); 
subplot(2,1,1); plot(tc, xc, 'LineWidth', 1); 
hold on;
stem(td1, xd1, '.'); 
hold off; grid on;
title('Sampling well above Nyquist (fs=8 kHz > 2*500 Hz)'); 


subplot(2,1,2); plot(tc, xc, 'LineWidth', 1); 
hold on;
stem(td2, xd2, '.'); 
hold off; grid on;
title('Sampling near Nyquist (fs=600 Hz ~ 2 * 500 Hz)'); 

%% Next -- 2.2 Deterministic vs Random
x_det = cos(2*pi*100*t);
x_rnd = filter([1 0.5],[1 -0.9], randn(size(t))); % colored random process

figure('Name','Deterministic vs Random');
subplot(2,1,1); plot(t(1:3000), x_det(1:3000)); 
title('Deterministic (100 Hz cosine)'); grid on;

subplot(2,1,2); plot(t(1:3000), x_rnd(1:3000)); 
title('Random (colored noise)'); grid on;

%% Next -- 2.3 Even/Odd Decomposition
x = sin(2*pi*5*t) + 0.3*cos(2*pi*7*t);
x_even = 0.5*(x + fliplr(x)); % fliplr ~ time reversal over finite grid
x_odd  = 0.5*(x - fliplr(x));

figure('Name','Even/Odd Decomposition');
subplot(3,1,1); plot(t, x); title('x(t)'); grid on;

subplot(3,1,2); plot(t, x_even); title('Even part'); grid on;

subplot(3,1,3); plot(t, x_odd);  title('Odd part');  grid on;

%% Next -- 2.4 Periodic vs Aperiodic
xP = cos(2*pi*25*t);                  % periodic 25 Hz cosine
xA = exp(-5*abs(t-0.5));              % more sharply localized aperiodic pulse

tview = t(1:5000);                    % view first 0.05 s for fair comparison

figure('Name','Periodic vs Aperiodic (Improved)');
subplot(2,1,1);
plot(tview, xP(1:length(tview)), 'LineWidth', 1.2);
title('Periodic: repeats indefinitely (25 Hz cosine)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(t, xA, 'LineWidth', 1.2);
title('Aperiodic: localized pulse (no repetition)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

%% Next -- 2.5 Energy vs Power Signals 

% Energy signal: 
xE = exp(-5*abs(t-0.5));      

% Power signal: 
xW = cos(2*pi*10*t);          

% Numerical approximations
E_xE = sum(abs(xE).^2)*Ts;    % approx ∫|xE|^2 dt
P_xW = mean(abs(xW).^2);      % approx limit_{T→∞} (1/T ∫|xW|^2 dt)

fprintf('Energy(xE) ≈ %.4f,   Avg Power(xW) ≈ %.4f\n', E_xE, P_xW);

% Visual comparison
tview = t(1:3000);    % view first 0.03 seconds

figure('Name','Energy vs Power Signals');
subplot(2,1,1);
plot(t, xE, 'LineWidth', 1.2);
title('Energy Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(tview, xW(1:length(tview)), 'LineWidth', 1.2);
title('Power Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

%%  -- 2.6 Time Shift & Time Scale
x0     = cos(2*pi*20*t);            % original 20 Hz
x_shift= cos(2*pi*20*(t-0.05));     % shift right by +0.05 s
x_scale= cos(2*pi*40*t);            % time-compressed (freq doubled to 40 Hz)

% Choose a window that actually shows shift (0.1 s is reasonable)
tview = t(t <= 0.1);
idx   = 1:length(tview);

figure('Name','Shift & Scale (Improved)');
subplot(3,1,1);
plot(tview, x0(idx), 'LineWidth', 1.2);
title('Original x(t) = cos(2\pi \cdot 20 t)');
xlabel('Time (s)'); grid on;

subplot(3,1,2);
plot(tview, x_shift(idx), 'LineWidth', 1.2);
title('Time-Shifted x(t-0.05): right shift by +0.05 s');
xlabel('Time (s)'); grid on;

subplot(3,1,3);
plot(tview, x_scale(idx), 'LineWidth', 1.2);
title('Time-Scaled (compressed): 40 Hz → twice as fast');
xlabel('Time (s)'); grid on;

%% 3) Classification of Systems (Ch. 2.4)
% Test of linearity; time-invariance (TI), Causality on simple examples.

%% 3.1 Linearity & Time Invariance 

% System S1: y[n] = 0.5*x[n] + x[n-1]
N  = 128; n = 0:N-1;
x1 = cos(2*pi*0.05*n);
x2 = sin(2*pi*0.10*n);
a = 1.2; b = -0.8;

S1 = @(x) 0.5*[0 x(1:end-1)] + x;   % simple FIR-like LTI system

%% --- Linearity Test ---
lhs = S1(a*x1 + b*x2);
rhs = a*S1(x1) + b*S1(x2);

figure('Name','Linearity Test for S1');
plot(n, lhs, 'LineWidth',1.2); hold on;
plot(n, rhs, '--', 'LineWidth',1.2);
legend('S1[a x_1 + b x_2]', 'a S1[x_1] + b S1[x_2]');
xlabel('n'); ylabel('Amplitude'); grid on;
title('Linearity Check (Curves Overlap if Linear)');

lin_err = norm(lhs-rhs);
fprintf('Linearity error (S1): ||lhs-rhs|| = %.3e\n', lin_err);

%% --- Time-Invariance Test ---
kShift = 3;
x_shiftTI = [zeros(1,kShift), x1(1:end-kShift)];
y1 = S1(x1);
y2 = S1(x_shiftTI);
y1_shift = [zeros(1,kShift), y1(1:end-kShift)];

figure('Name','Time-Invariance Test for S1');
plot(n, y2, 'LineWidth',1.2); hold on;
plot(n, y1_shift, '--', 'LineWidth',1.2);
legend('S1[x(n-k)]', 'S1[x(n)] shifted by k');
xlabel('n'); ylabel('Amplitude'); grid on;
title(sprintf('Time Shift Check (k=%d)', kShift));

ti_err = norm(y2 - y1_shift);
fprintf('Time-invariance error (S1): ||y2 - y1_shift|| = %.3e\n', ti_err);


%% 3.3 Convolution Demo (Impulse Response & Frequency Response)
% Goal:
%   Show that a lowpass FIR filter passes low-frequency components
%   (20 Hz) and attenuates higher-frequency components (60 Hz).

% Input signal: sum of two cosines
f1 = 20;      % low frequency (Hz)
f2 = 60;      % higher frequency (Hz)
x_in = cos(2*pi*f1*t) + 0.5*cos(2*pi*f2*t);

% Design a simple lowpass FIR via windowed-sinc
% Place cutoff between 20 Hz and 60 Hz, say at 40 Hz
fc   = 40;             % cutoff frequency in Hz
fcN  = fc / (fs/2);    % normalized to Nyquist (0..1)
Nf   = 201;            % filter length (odd for symmetry)
m    = -(Nf-1)/2:(Nf-1)/2;
h_id = 2*fcN*sinc(2*fcN*m);     % ideal LPF impulse response
h    = h_id .* hann(Nf)';       % apply Hann window

% Filter the input
y_out = conv(x_in, h, 'same');

% Time-domain view (short window for clarity)
idx = 1:4000;           % show first 0.04 s
tview = t(idx);

figure('Name','Convolution Demo: LPF on 20+60 Hz');
subplot(3,1,1);
plot(tview, x_in(idx), 'LineWidth', 1.1);
title('Input: 20 Hz + 0.5·60 Hz'); grid on;
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,2);
plot(tview, y_out(idx), 'LineWidth', 1.1);
title('Output after Lowpass Filter (Cutoff ≈ 40 Hz)'); grid on;
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,3);
stem(m, h, 'filled'); grid on;
title('Impulse Response h[n]'); xlabel('n'); ylabel('h[n]');

% Frequency response and tone markers
figure('Name','Frequency Response of LPF');
[Hf, w] = freqz(h, 1, 2048, fs);   % freq response in Hz
plot(w, 20*log10(abs(Hf)+1e-12), 'LineWidth',1.1); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('LPF Magnitude Response');

hold on;
yL = ylim;
plot([f1 f1], yL, 'r--', 'LineWidth',1);   % mark 20 Hz
plot([f2 f2], yL, 'g--', 'LineWidth',1);   % mark 60 Hz
legend('LPF |H(f)|', '20 Hz tone', '60 Hz tone', 'Location','Best');
hold off;

%% 4) Signal Representation (SDR4Engineers §2.3): Complex Baseband & I/Q
% Build I/Q baseband, upconvert to RF, then downconvert back to baseband.
% Also fix plotting so traces are clearly visible (no "solid stains").

%% Next -- 4.1 I/Q Baseband and Upconversion
fbb   = 1e3;                             % baseband tone frequency (1 kHz)
Ibb   = cos(2*pi*fbb*t);
Qbb   = 0.7*sin(2*pi*fbb*t + pi/6);
ybb   = Ibb + 1j*Qbb;
fc_rf = 20e3;                            % RF carrier 20 kHz
wc    = 2*pi*fc_rf;
r     = Ibb.*cos(wc*t) - Qbb.*sin(wc*t);

% Visualize spectra (optional quick look with Welch PSD)
figure('Name','PSD Views (Optional)');
pwelch(ybb, [], [], [], fs, 'centered');
title('Baseband Complex Signal PSD');

%% Next -- 4.2 Downconversion & Proper Plotting (Fixed)
% Mix down and lowpass to recover I and Q.
% Choose a sensible LPF cutoff that passes fbb=1 kHz but rejects 2*fc terms.
Fcut  = 5e3;                  % 5 kHz cutoff
Ir    = 2*lowpass(r.*cos(wc*t), Fcut, fs);   % factor 2 compensates 1/2 from product
Qr    = -2*lowpass(r.*sin(wc*t), Fcut, fs);

% For CLEAR plots: show a short segment and use different styles
idx = 1:2000;                  % show first 2000 samples (~0.02 s)
tseg = t(idx);
Iseg_true = Ibb(idx);  Iseg_rec  = Ir(idx);
Qseg_true = Qbb(idx);  Qseg_rec  = Qr(idx);

figure('Name','I Component (Segment)');
plot(tseg, Iseg_true, 'LineWidth', 1.2); hold on;
plot(tseg, Iseg_rec,  '--', 'LineWidth', 1.2);
grid on; legend('I true','I rec'); 
xlabel('Time (s)'); ylabel('Amplitude'); 
title('I(t): true vs recovered'); hold off;

figure('Name','Q Component (Segment)');
plot(tseg, Qseg_true, 'LineWidth', 1.2); hold on;
plot(tseg, Qseg_rec,  '--', 'LineWidth', 1.2);
grid on; legend('Q true','Q rec'); 
xlabel('Time (s)'); ylabel('Amplitude'); 
title('Q(t): true vs recovered'); hold off;

% Quantitative agreement
rmseI = sqrt(mean((Ibb - Ir).^2));
rmseQ = sqrt(mean((Qbb - Qr).^2));
fprintf('RMSE(I): %.3e,  RMSE(Q): %.3e\n', rmseI, rmseQ);

%% Next -- 4.3 LO Phase Error Demonstration
phi_err = deg2rad(5); % 5 degrees phase error at RX
Ir_pe   = 2*lowpass(r.*cos(wc*t + phi_err), Fcut, fs);
Qr_pe   = -2*lowpass(r.*sin(wc*t + phi_err), Fcut, fs);
idx2 = 1:2000; tseg2 = t(idx2);

figure('Name','LO Phase Error Effect');
plot(tseg2, Ibb(idx2), 'LineWidth', 1.2); hold on;
plot(tseg2, Ir_pe(idx2), '--', 'LineWidth', 1.2);
grid on; 
legend('I true','I rec (5° error)'); 
title('I-path: true vs recovered with 5 deg LO phase error'); hold off;

% Q16: What happens in the constellation view with LO phase error?

%% 5) Information Theory (SDR4Engineers Ch. 1; Principles §1.5): Capacity Intuition
% Shannon capacity: C = B * log2(1 + SNR).

%% 5.1 Capacity vs SNR for Fixed B
B = 20e3;                                       % 20 kHz
SNRdB = -10:2:30; 
SNRlin = 10.^(SNRdB/10);
C = B*log2(1+SNRlin);

figure('Name','Capacity vs SNR');
plot(SNRdB, C/1e3, 'o-'); grid on;
xlabel('SNR (dB)'); 
ylabel('Capacity (kb/s)'); 
title(sprintf('C vs SNR (B = %.0f kHz)', B/1e3));

% Q17: At very low SNR, which is more effective: +3 dB SNR or doubling B?

%% Next -- 5.2 Capacity Heatmap: Tradeoff of B and SNR
Bvals   = linspace(5e3, 200e3, 20);
SNRdBs  = linspace(-5, 25, 19);
[BB, SS] = meshgrid(Bvals, 10.^(SNRdBs/10));
CH = BB.*log2(1+SS)/1e3; % kb/s

figure('Name','Capacity Heatmap');
imagesc(Bvals/1e3, SNRdBs, CH); axis xy;
xlabel('Bandwidth (kHz)'); ylabel('SNR (dB)'); colorbar; 
title('Capacity (kb/s)');

% Q18: Why can capacity scale almost linearly with B at low SNR?

%% 6) Mini Exercises (with Instructor Notes)
% (a) Change fc_rf to 200 kHz and fs to 1 MHz; verify recovery (expect similar RMSE).
% (b) Replace Ibb/Qbb with RRC-shaped BPSK symbols; add 3° LO phase error and measure SER.
% (c) For fixed SNR=10 dB, plot capacity vs bandwidth 5–200 kHz; identify diminishing returns when system cannot exploit more B.

%% 7) References Mapping (for your prep)
% Principles of Modern Communication Systems
%  - 1.1–1.5: system model, frequency allocation, intro to information theory
%  - 2.1–2.4: signal and system classifications, LTI fundamentals
% SDR for Engineers
%  - 1: capacity intuition in context of SDR workflows
%  - 2.3: complex baseband (I/Q), frequency conversion, and IQ impairments

%% Local Function Utilities (used in some sections if you extend)
% (Add any helper functions below; Live Editor supports local functions at end of script.)
