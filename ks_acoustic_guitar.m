% -------------------------------------------------------------------------
% Script that simulates acoustic guitar sound using modified Karplus-Strong 
% algorithm.
%
% This script includes following extensions to Karplus-Strong algorithm:
% 1. Decay stretching [1].
% 2. User-specified T60 time for loss parameter calculation [1]. T60 and 
%    note duration are controlled independently to allow for more timbre
%    variations (especially for synthesising melodies where notes can be 
%    muted).
% 3. Fix for the discontinuity in the end of notes [1].
% 4. Different excitation techniques such as pluck, slide, bend, hammer-on,
%    pull-off, vibrato [1].
% 5. Pick position simulation [1].
% 6. Pick "up" and "down" directions simulation [1].
% 7. Dual-polarity string model [2].
% 8. Convolution with a guitar body impulse response [1,3]. The convolution
%    is performed only if the sampling rate of the algorithm and the 
%    impulse response are equal.
%
% A simple sequencer is used to synthesise a guitar melody and demonstrate
% the developed physical model.
%
% NB: body impulse response is not included in the repository and can be
% downloaded from [3].
%
% REFERENCES:
% [1] Jaffe, D. A., & Smith, J. O. (1983). Extensions of the Karplus-Strong
%     Plucked-String Algorithm. Computer Music Journal, 7(2), 56-69. 
%     https://doi.org/10.2307/3680063
% [2] Karjalainen, M., Valimaki, V., & Tolonen, T. (1998). Plucked-String 
%     Models: From the Karplus-Strong Algorithm to Digital Waveguides and 
%     beyond. Computer Music Journal, 22(3), 17-32. 
%     https://doi.org/10.2307/3681155
% [3] Acoustic guitars impulse responses database.
%     URL (general): http://acousticir.free.fr
%     URL (Taylor_220ce_K_DLX): https://drive.google.com/drive/u/0/folders/1OznSBjCeBcWhK77hdS5U7LUFa8Ml-l06
%
% Author: Victor Zheleznov
% Date: 31/01/2023
% -------------------------------------------------------------------------

% clear workspace
clear all; close all;

% define parameters
Fs = 48000; % sample rate [Hz]
bpm = 100;  % sequencer bpm
st = -5;    % sequencer pitch transposition [semitones]
p = 0.9;    % max coef. value for lowpass filter for pick direction simulation

% parameters for slightly mistuned model (for polarisation simulation)
par_pol.dt60 = 5; % max difference in T60 time [%]
par_pol.df0 = 4;  % max difference in fundamental frequency [cents]
par_pol.w = 0.2;  % weight for mistuned model

% parameters for the end of the note
par_end.t_end = 5e-3;  % when to change the loss factor [sec] (counting from the end)
par_end.factor = 1e-1; % by which factor to multiply the loss factor

% load guitar body impulse response
[h, hFs] = audioread("Taylor_220ce_K_DLX_48000.wav"); % can be downloaded from [3]

% define melody
notes = struct('f0',   {247*2^(10/12), 330*2^(10/12), 330*2^(8/12), 247*2^(10/12), 330*2^(8/12)},... % note pitch [Hz] (f0 = 0 mean a rest)
               'dur',  {1/12,          1/6,           1/12,         1/6,           1           },... % note relative duration
               'mode', {"pluck",       "slide",       "pluck",      "vibrato",     "bend"      },... % note excitation type: pluck, bend, slide, hammer, vibrato
               'df',   {0,             200,           0,            30,            100         },... % frequency shift for a note [cents] (can be negative i.e. for "hammer" mode to simulate pull-off) 
               't_df', {0,             0.2,           0,            0.95,          0.4         },... % frequency shift duration relative to the note length (for "hammer" mode specifies when to trigger it)
               'R',    {0.9,           0.5,           0.6,          0.4,           0.6         },... % dynamics coefficient in (0,1) (from high to low dynamics)
               'S',    {0.5,           0.6,           0.5,          0.6,           0.8         },... % decay stretching factor in (0,1) (minimum decay for S = 0.5)
               'nu',   {0.3,           0.3,           0.4,          0.3,           0.4         },... % pick position in (0,1) (fraction of the string between the bridge and pluck point)
               't60',  {1.5,           2.5,           1.5,          2.5,           2.4         },... % desired t60 time [sec]
               'fv',   {0,             0,             0,            6,             0           });   % vibrato frequency [Hz]

% generate alternating picking pattern
p_pattern = num2cell(p*mod(0:length(notes)-1,2));
[notes.p] = p_pattern{:};

% convert note duration to seconds
dur_sec = num2cell([notes.dur]*4*60/bpm);
[notes.dur] = dur_sec{:};

% transpose notes
f0_tr = num2cell([notes.f0]*2^(st/12));
[notes.f0] = f0_tr{:};

% synthesise melody
M = floor(Fs*[dur_sec{:}]);
y = zeros(sum(M),1);
idx = 1;
for i = 1:length(notes)
    y(idx:idx+M(i)-1) = synthesise_acoustic_guitar_tone(notes(i), par_pol, par_end, Fs);
    idx = idx + M(i);
end

% convolve output with guitar body impulse response and normalise
if hFs == Fs
    y = myfastconv(y,h);
else
    y = y./max(abs(y));
end

% listen to the output, save audio and plot spectogram
soundsc(y, Fs);
audiowrite('ks_acoustic_guitar.wav', y, Fs);
myspec(y, Fs, 2048, 0.75);

%% FUNCTIONS
% synthesise acoustic guitar tone
% input:
%   note - struct with note information;
%   par_pol - struct with parameters for second polarisation simulation; 
%   par_end - struct with parameters for the end of the note;
%   Fs - sampling rate [Hz].
% output:
%   y - synthesised waveform.
function y = synthesise_acoustic_guitar_tone(note, par_pol, par_end, Fs)
    % save data to local variables
    f0 = note.f0;         % fundamental frequency [Hz]
    dur = note.dur;       % note duration [sec]
    mode = note.mode;     % note excitation type: pluck, bend, slide, hammer, vibrato
    df = note.df;         % frequency shift for a note [cents]
    t_df = note.t_df*dur; % frequency shift duration [sec]
    R = note.R;           % dynamics coefficient
    S = note.S;           % decay stretching factor 
    p = note.p;           % "up" or "down" picking
    nu = note.nu;         % pick position
    t60 = note.t60;       % T60 [sec]
    fv = note.fv;         % vibrato frequency [Hz]
    df0 = par_pol.df0;    % max difference in fundamental frequency [cents]
    dt60 = par_pol.dt60;  % max difference in t60 [%]
    w = par_pol.w;        % weight for mistuned model
    
    % check if input is a rest
    M = floor(Fs*dur); % note duration in samples
    if f0 == 0
        y = zeros(M,1);
        return
    end
    
    % generate mistuned parameters (for polarisation simulation)
    f0_ = f0*2^((-df0 + 2*df0*rand)/1200);
    t60_ = t60*(1 + 1e-2*dt60*rand*sign(f0-f0_)); 
    
    % generate time-dependant note frequency
    Md = floor(Fs*t_df);                                                               % frequency modulation duration in samples
    N0 = max(floor(Fs/f0 - S), floor(Fs/f0_ - S));                                     % max initial delay line size
    Nmax = max([N0, floor(Fs/(f0*2^(df/1200)) - S), floor(Fs/(f0_*2^(df/1200)) - S)]); % max overall delay line size
    if mode == "pluck"
		f = f0*ones(M,1);
	elseif mode == "bend"
		c = [zeros(Nmax+1,1); (df/Md:df/Md:df).'; df*ones(M-Md-Nmax-1,1)];
		f = f0 * 2.^(c/1200);
	elseif mode == "vibrato"
		arg = 2*pi*fv*(0:1/Fs:floor(t_df*fv)/fv).';
		L = length(arg);
		win = 0.5*(1 - cos(2*pi.*(0:L-1)./L)).';
		c = [zeros(Nmax+1,1); df*win.*sin(arg); zeros(M-L-Nmax-1,1)];
		f = f0 * 2.^(c/1200);
	elseif mode == "hammer"
        L = max(Nmax+1, Md);
		c = df.*[zeros(L, 1); ones(M-L, 1)];
		f = f0 * 2.^(c/1200);
	elseif mode == "slide"
		c = [zeros(Nmax+1,1); (df/Md:df/Md:df).'; df*ones(M-Md-Nmax-1,1)];
		c = floor(c./100);
		f = f0 * 2.^(c/12);
    end
	
    % excitation model
    exc_signal = excitation_model(N0, R, nu, p);
    
    % string model (two polarisations)
    f = [f, f.*(f0_/f0)];
    t60 = [t60, t60_];
    w = [1-w, w];
    y = zeros(M,1);
    for i = 1:2
        % calculate derived parameters
        Nexact = Fs./f(:,i) - S;                             % delay line length
        N = floor(Nexact);                                   % truncated delay line length
        P = Nexact - N;                                      % fractional delay
        C = (1-P)./(1+P);                                    % tuning all-pass filter coefficient
        tau = t60(i)/log(1000);
        rho = (1/abs(cos(pi*f(1,i)/Fs)))*exp(-1/f(1,i)/tau); % loss coefficient
        % synthesise string using extended Karplus–Strong algorithm
        y = y + w(i)*acoustic_guitar_string_model(exc_signal, M, N, C, S, rho, par_end, Fs);
    end
end

% create excitation signal using noise generator
% input:
%   N0 - desired number of samples (output length is N0+1);
%   R - dynamics coefficient in (0,1);
%   nu - relative pick position;
%   p - pick direction coefficient in [0,1).
% output:
%   yf - synthesised excitation signal.
function yf = excitation_model(N0, R, nu, p)
    % generate noise
	v = -1 + 2*rand(N0+1,1); % input for dynamics filter
	
	% apply lowpass filter for dynamics simulation
	yd = filter(1-R, [1 -R], v);
	
	% apply comb filter for pick position simulation
	ye = filter([1, zeros(1, floor(nu*N0)-1), -1], 1, yd);
	
	% apply lowpass filter for pick direction simulation
	yf = filter(1-p, [1 -p], ye);
end

% simulate acoustic guitar string model
% input:
%   exc_signal - excitation signal;
%   M - output vector length [samples];
%   N - truncated delay line length;
%   C - fractional delay;
%   S - decay stretching factor in (0,1)
%   rho - loss factor;
%   par_end - struct with parameters for the end of the note;
%   Fs - sampling frequency [Hz].
% output:
%   y - synthesised string output.
function y = acoustic_guitar_string_model(exc_signal, M, N, C, S, rho, par_end, Fs)
    % parameters for the end of the note
    t_end = par_end.t_end;   % when to change the loss factor [sec]
    factor = par_end.factor; % by which factor to multiply the loss factor
    
    % initialise output array
    y = zeros(M,1);          
    y(1:N(1)+1) = exc_signal(1:N(1)+1);
	
    % process the Karplus-Strong algorithm
	yp1 = 0;                                           % delayed output from the tuning filter
    for n = N(1)+1:M-1
		yp0 = C(n)*y(n-N(n)+1) + y(n-N(n)) - C(n)*yp1; % allpass filter for tuning 
		y(n+1) = rho*((1-S)*yp0 + S*yp1);              % lowpass filter for damping
		yp1 = yp0;                                     % save data to the delayed sample
        % decrease loss value at the end of the note
        if n == M-1-round(t_end*Fs)
            rho = rho*factor;
        end
    end
end

% perform convolution of two real (!) vectors utilising their fft()
% input:
%   h - first vector;
%   x - second vector;
% ouput:
%   y - normalised (in [-1,1]) convolution result, 
%       length(y) = length(h) + length(x) - 1.
function y = myfastconv(h,x)
    % guarantee column vectors
    h = h(:);
    x = x(:);
    
    % determine fft size
    N = length(h) + length(x) - 1;
    NFFT = 2^(ceil(log(N)/log(2))); % next power of 2
    NFFT_2 = NFFT / 2 + 1;
    
    % calculate convolution
    XF = fft(x, NFFT);
    XF = XF(1:NFFT_2);
    HF = fft(h, NFFT);
    HF = HF(1:NFFT_2);
    YF = HF.*XF;
    
    % make DC and Nyquist bins real because due to the machine floating
    % point error they can have small imaginary parts (because sin(N*pi) 
    % will not be zero where N is an integer)
    YF(1) = real(YF(1));
    YF(end) = real(YF(end));
    
    % take inverse transform
    YF = [YF; conj(YF(end-1:-1:2))];
    y = ifft(YF);
    y = y(1:N);
    
    % normalise output
    y = y./max(abs(y));
end

% create a spectogram plot of an input signal
% input:
%   x - mono input signal;
%   Fs - sampling frequency [Hz];
%   N - frame length;
%   O - overlap factor (between 0 and 1).
function [] = myspec(x, Fs, N, O)
    % find hop size
    HA = round(N - O*N);

    % generate window
    win = 0.5*(1 - cos(2*pi.*(0:N-1)./N)).';

    % calculate number of frames
    L = length(x);
    NF = ceil(L/HA);
    x = [x; zeros((NF-1)*HA+N-L,1)];
    
    % STFT size
    NFFT = 2^(ceil(log(N)/log(2))); % next power of 2
    NFFT_2 = NFFT / 2 + 1;

    % calculate STFT
    STFT = zeros(NFFT_2, NF);
    for m = 0:NF-1
        x_frame = win.*x((1:N).'+m*HA);
        X = fft(x_frame, NFFT);
        STFT(:,m+1) = X(1:NFFT_2);
    end
    
    % plot spectogram
    fig_spec = figure;
    t = ((0:NF-1).*HA/Fs).';
    freq = (0:Fs/NFFT:Fs/2).';
    STFT_dB = 20*log10(abs(STFT));
    max_dB = max(max(STFT_dB));
    imagesc(t, freq, STFT_dB, 'CDataMapping', 'scaled');
    c = colorbar;
    c.Label.String = 'dB';
    colormap hot
    caxis([max_dB-60, max_dB]);
    xlim([0 t(end)]);
    ylim([0 freq(end)]);
    ax_spec = fig_spec.CurrentAxes;
    set(ax_spec, 'YDir', 'normal');
    set(ax_spec, 'YTick', 0:1000:Fs/2);
    set(ax_spec, 'YTickLabel', 0:1000:Fs/2);
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('Frequency [Hz]', 'interpreter', 'latex');
    title_str = sprintf("Spectogram with frame length = $%d$ ms and overlap factor = $%d$\\%%", floor((N/Fs)*1e3), O*1e2);
    title(title_str, 'interpreter', 'latex');
end