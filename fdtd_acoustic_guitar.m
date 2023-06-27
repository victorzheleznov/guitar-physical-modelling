% -------------------------------------------------------------------------
% Script for simulating an acoustic guitar using finite-difference
% time-domain method. 
%
% Simulation includes the following extensions of the stiff string model:
% 1. Frequency dependent loss [2].
% 2. Tension modulation effect [1].
% 3. Fretboard collision model [1].
% 4. Finger collision model [1].
% 5. Six guitar strings with realistic physical parameters and adjustable
%    tuning.
% 6. Randomisation of excitation and finger parameters.
% 7. String velocity output option.
% 8. Stereo signal emulation.
%
% Parameters for guitar strings were obtained from different sources:
% 1. Strings length and radius: 
%    https://www.elixirstrings.com/support/string-tension-for-tuning-guitar
%    (Acoustic Phosphor Bronze with NANOWEB Coating, Medium Strings)
% 2. Phosphor Bronze density and Young's modulus:
%    https://www.azom.com/article.aspx?ArticleID=6415
% 3. T60 at different frequencies were obtained using my personal guitar 
%    recordings and Room EQ Wizard software:
%    https://www.roomeqwizard.com/
%
% NB: since simulation time is long due to the nonlinearities in the model,
% the demo recordings are included in the repository:
% 1. fdtd_acoustic_guitar_open_strings.wav
%    Six open strings are exited by a pluck of increasing amplitude (from 
%    1 N for low E string to 3.5 N for high E string). In the end you can
%    clearly hear string collision with the fretboard.
% 2. fdtd_acoustic_guitar_chord_tapping.wav
%    Six string are tapped against the fretboard by a finger to produce 
%    a Fmaj7 chord. You can hear rattling from finger/string/fretboard 
%    collisions.
%
% By default the program runs with collisions and plotting turned on to
% demonstrate the developed model by tapping strings against fretboard. 
% The specified configuration corresponds to the demo recording
% fdtd_acoustic_guitar_chord_tapping.wav.
%
% You can turn collisions and plotting off (by using 
% flags.calc_collisions = 0 and flags.plot_on = 0) to excecute the program 
% much faster and hear the plucking of open strings. You'll need to specify 
% 'feamp' parameter in the 'par_io' struct to set the excitation amplitude 
% not to zero.
%
% References:
% [1] S.Bilbao, and  A.Torin, "Numerical Modeling and Sound Synthesis for 
%     Articulated String/Fretboard Interactions," J. Audio Eng. Soc., vol. 
%     63, no. 5, pp. 336-347, (2015 May.).
%     https://doi.org/10.17743/jaes.2015.0023
% [2] Stefan Bilbao. Numerical Sound Synthesis: Finite Difference Schemes 
%     and Simulation in Musical Acoustics.
%     https://doi.org/10.1002/9780470749012
% 
% Author: Victor Zheleznov
% Date: 13/02/2023
% -------------------------------------------------------------------------

clc; clear all; close all

% flags
flags.plot_on = 1;         % in-loop plotting on (1) or off (0)
flags.itype = 1;           % type of input: 1 - pluck, 2 - strike
flags.otype = 2;           % output type: 1 - displacement, 2 - velocity
flags.calc_collisions = 1; % calculate collisions with fretboard and finger: on (1) or off (0)

% physical strings parameters
%                            high E     A          D          G          B          low E
par_phys = struct('r',      {0.0001651, 0.0002159, 0.0003302, 0.0004445, 0.0005715, 0.0007112},... % radius [m]
                  'rho',    {8860,      8860,      8860,      8860,      8860,      8860     },... % density [kg/m^3]
                  'L',      {0.6477,    0.6477,    0.6477,    0.6477,    0.6477,    0.6477   },... % length [m]
                  'E',      {117e9,     117e9,     117e9,     117e9,     117e9,     117e9    },... % Young's modulus [Pa]
                  'f0',     {329.63,    246.94,    196.00,    146.83,    110.00,    82.41    },... % fundamental frequency [Hz] (currently in EADGBE tuning)
                  'T60_f1', {13,        13,        16.9,      14.8,      24.7,      17.8     },... % T60 [sec] at frequency f1
                  'f1',     {333.7,     250,       198.4,     148.7,     111.4,     81.05    },... % frequency f1 [Hz]
                  'T60_f2', {8,         7.9,       14.5,      13.5,      13.2,      13.9     },... % T60 [sec] at frequency f2
                  'f2',     {667.4,     749.2,     369.9,     432.8,     216.4,     250      });   % frequency f2 [Hz]

% input/output
SR = 44100;
Tf = 2;
par_io = struct('SR',      {SR                                      },... % sample rate [Hz]
                'Tf',      {Tf                                      },... % duration of simulation [sec]
                'xi',      {0.9                                     },... % coordinate of excitation (normalised, 0-1)
                'feamp',   {0                                       },... % peak amplitude of excitation [N]
                'exc_dur', {0.001                                   },... % duration of excitation [sec]
                'exc_st',  {0.1,   0.07,  0.05,  0.03,  0.01,  0.001},... % start time of excitation [sec]
                'xo',      {0.85                                    },... % coordinate of output (normalised, 0-1)
                'nf',      {5,     5,     5,     7,     8,     5    },... % fret number (for finger position)
                'f0amp',   {10,    10,    10,    15,    15,    15   },... % finger force amplitude [N]
                'fin_dur', {20                                      },... % finger force duration [sec]
                'fin_st',  {0.495, 0.395, 0.295, 0.195, 0.095, 0    });   % finger force start time [sec]

% fretboard
par_fret.Nfret = 20;    % number of frets
par_fret.Hfret = 5e-4;  % frets height relative to the fretboard [m]
par_fret.Hback = -1e-3; % fretboard position relative to the string [m]
par_fret.Kc = 1e15;     % scaling parameter for fretboard collision potential
par_fret.ac = 2.3;      % degree parameter for fretboard collision potential

% finger
par_fing.Mf = 5e-3;  % finger mass [kg]
par_fing.Kf = 1e10;  % scaling parameter for finger collision potential
par_fing.af = 2.3;   % degree parameter for finger collision potential
par_fing.bf = 1;     % loss parameter for finger collision potential
par_fing.uf0 = 2e-3; % finger initial position [m]
par_fing.vf0 = 0;    % finger initial velocity [m]

% Newton Raphson
par_nr.tol_nr = 1e-8; % tolerance value for Newton-Raphson algorithm
par_nr.max_iter = 20; % max iterations for Newton-Raphson algorithm

% randomisation
par_rand.exc_st_perc = 1e-1;   % percentage for excitation start randomisation
par_rand.xi_perc = 5e-2;       % percentage for excitation position randomisation
par_rand.feamp_perc = 1e-1;    % percentage for excitation force amplitude randomisation
par_rand.fin_st_perc = 5e-2;   % percentage for finger start randomisation
par_rand.xf_perc = 5e-2;       % percentage for finger position randomisation
par_rand.f0amp_perc = 1e-1;    % percentage for finger force amplitude randomisation

% output parameters
wl = [0.35; 0.4; 0.48; 0.54; 0.59; 0.65];
wr = 1-wl;
fade_perc = 1e-1;

% simulate strings
Nf = floor(SR*Tf);
y = zeros(Nf,6);
for i = 1:6
    disp("string " + i + " simulation...");
    y(:,i) = string_fdtd(par_phys(i), par_io(i), par_fret, par_fing, par_nr, par_rand, flags);
end

% create stereo signal
y_stereo = zeros(Nf,2);
Nfade = floor(fade_perc*Nf); % fade out
win = [ones(Nf-Nfade,1); (Nfade-1:-1:0).'/Nfade];
y_stereo(:,1) = win.*(y*wl);
y_stereo(:,2) = win.*(y*wr);

% listen
soundsc(y_stereo, SR);

% save
audiowrite('fdtd_acoustic_guitar.wav', y_stereo/max(abs(y_stereo), [], 'all'), SR);
save('fdtd_acoustic_guitar.mat');

% plot spectogram
myspec(y_stereo(:,1), SR, 2048, 0.75);

%% FUNCTIONS
function y = string_fdtd(par_phys, par_io, par_fret, par_fing, par_nr, par_rand, flags)
    % flags
    plot_on = flags.plot_on; % in-loop plotting on (1) or off (0)
    itype = flags.itype;     % type of input: 1 - pluck, 2 - strike
    otype = flags.otype;     % output type: 1 - displacement, 2 - velocity
    calc_collisions = flags.calc_collisions;

    % physical parameters
    f0 = par_phys.f0;         % fundamental frequency [Hz]
    r = par_phys.r;           % string radius [m]
    E = par_phys.E;           % Young's modulus [Pa]
    rho = par_phys.rho;       % density [kg/m^3]
    T60_f1 = par_phys.T60_f1; % T60 [sec] at frequency f1
    f1 = par_phys.f1;         % frequency f1 [Hz]
    T60_f2 = par_phys.T60_f2; % T60 [sec] at frequency f2
    f2 = par_phys.f2;         % frequency f2 [Hz]
    L = par_phys.L;           % length [m]
    
    % input/output parametrs
    SR = par_io.SR;           % sample rate [Hz]
    Tf = par_io.Tf;           % duration of simulation [sec]
    xi = par_io.xi;           % coordinate of excitation (normalised, 0-1)
    feamp = par_io.feamp;     % peak amplitude of excitation [N]
    exc_dur = par_io.exc_dur; % duration of excitation [sec]
    exc_st = par_io.exc_st;   % start time of excitation [sec]
    xo = par_io.xo;           % coordinate of output (normalised, 0-1)
    nf = par_io.nf;           % finger fret number
    f0amp = par_io.f0amp;     % finger force amplitude [N]
    fin_dur = par_io.fin_dur; % duration of excitation [sec]
    fin_st = par_io.fin_st;   % start time of excitation [sec]
    
    % fretboard parameters
    Nfret = par_fret.Nfret; % number of frets
    Hfret = par_fret.Hfret; % frets height relative to the fretboard [m]
    Hback = par_fret.Hback; % fretboard position relative to the string [m]
    Kc = par_fret.Kc;       % scaling parameter for fretboard collision potential
    ac = par_fret.ac;       % degree parameter for fretboard collision potential
    
    % finger
    Mf = par_fing.Mf;   % finger mass [kg]
    Kf = par_fing.Kf;   % scaling parameter for finger collision potential
    af = par_fing.af;   % degree parameter for finger collision potential
    bf = par_fing.bf;   % loss parameter for finger collision potential
    uf0 = par_fing.uf0; % finger initial position [m]
    vf0 = par_fing.vf0; % finger initial velocity [m]
    
    % Newton Raphson
    tol_nr = par_nr.tol_nr;     % tolerance value for Newton-Raphson algorithm
    max_iter = par_nr.max_iter; % max iterations for Newton-Raphson algorithm
    
    % randomisation
    exc_st_perc = par_rand.exc_st_perc; % percentage for excitation start randomisation
    xi_perc = par_rand.xi_perc;         % percentage for excitation position randomisation
    feamp_perc = par_rand.feamp_perc;   % percentage for excitation force amplitude randomisation
    fin_st_perc = par_rand.fin_st_perc; % percentage for finger start randomisation
    xf_perc = par_rand.xf_perc;         % percentage for finger position randomisation
    f0amp_perc = par_rand.f0amp_perc;   % percentage for finder force amplitude randomisation
    
    % check parameters
    assert(itype == 1 || itype == 2, 'Type of input should be 1 or 2!');
    assert(otype == 1 || otype == 2, 'Type of output should be 1 or 2!');
    assert(r > 0, 'String radius should be positive!');
    assert(E >= 0, "Young's modulus should be non-negative!");
    assert(rho > 0, 'Density should be positive!');
    assert(T60_f1 > 0 && T60_f2 > 0, 'T60 should be positive!');
    assert((f1 > 0 && f2 > 0) && (f1 < SR/2 && f2 < SR/2), 'Frequencies for T60 should be within Nyquist range!');
    assert(L > 0, 'Length should be positive!');
    assert(SR > 0, 'Sample rate should be positive!');
    assert(Tf > 0, 'Duration of simulation should be positive!');
    assert((xi >= 0) && (xi <= 1), 'Normalised coordinate of excitation should be in [0,1]!');
    assert((xo >= 0) && (xo <= 1), 'Normalised coordinate of output should be in [0,1]!');
    assert(exc_dur > 0, 'Excitation duration should be positive!');
    assert(exc_st > 0, 'Start time of excitation duration should be positive!');
    assert(exc_st + exc_dur < Tf, 'Excitation duration should be shorter than simulation duration!');
    assert((Nfret > 0) && (mod(Nfret,1) == 0), 'Number of frets should be positive and integer!');
    assert((Kc > 0) && (Kf > 0), 'Scaling parameters for collision potentials should be positive!');
    assert((ac > 1) && (af > 1), 'Degree parameters for collision potentials should larger than 1!');
    assert(Mf > 0, 'Finger mass should be positive!');
    assert(bf > 0, 'Loss parameter for finger collision potential should be positive!');
    assert((nf > 0) && (nf < Nfret), 'Finger fret number should be positive and less than number of frets!');
    
    % derived parameters
    A = pi*r^2;            % string cross-sectional area
    T = 4*L^2*f0^2*rho*A;  % string tension
    I = 0.25*pi*r^4;       % string moment of inertia
    c = sqrt(T/(rho*A));   % wave speed
    K = sqrt(E*I/(rho*A)); % stiffness constant 
    k = 1/SR;              % time step
    Nf = floor(SR*Tf);     % number of time steps
    
    % loss
    xi_f1 = xi_f(f1, c, K);
    xi_f2 = xi_f(f2, c, K);
    sig0 = 6*log(10)/(xi_f2-xi_f1)*(xi_f2/T60_f1 - xi_f1/T60_f2); % frequency independent loss parameter
    sig1 = 6*log(10)/(xi_f2-xi_f1)*(-1/T60_f1 + 1/T60_f2);        % frequency dependent loss parameter
    
    % calculate grid spacing
    hmin = sqrt(0.5*k*(c^2*k + 4*sig1 + sqrt((c^2*k+4*sig1)^2 + 16*K^2)));
    N = floor(L/hmin);
    h = L/N;
    
    % check end points
    assert(N < 10000, 'Number of segments should be less than 10000!');
    assert((xo*L >= h) && (xo*L <= L-h), 'Coordinate of output should be at least h metres away from endpoints!');
    assert((xi*L >= h) && (xi*L <= L-h), 'Coordinate of excitation should be at least h metres away from endpoints!');
    
    % calculate pluck or strike
    exc_st = exc_st*(1 + exc_st_perc*(-1 + 2*rand));
    feamp = feamp*(1 + feamp_perc*(-1 + 2*rand));
    n_exc_st = floor(exc_st/k);
    n_exc_dur = floor(exc_dur/k);
    fe = zeros(Nf,1);
    fe(n_exc_st:n_exc_st + n_exc_dur) = 0.5*feamp*(1 - cos(itype*pi*(0:n_exc_dur)/n_exc_dur));
    
    % calculate matrix form
    e = ones(N-1,1);
    Dxx = spdiags([e -2*e e], -1:1, N-1,N-1)./h^2;
    Dxxxx = spdiags([e -4*e 6*e -4*e e], -2:2, N-1,N-1);
    Dxxxx(1,1) = 5;
    Dxxxx(N-1,N-1) = 5;
    Dxxxx = Dxxxx./h^4;
    B = 2*speye(N-1) + (k^2*c^2 + 2*sig1*k)*Dxx - k^2*K^2*Dxxxx;
    
    % calculate excitation position
    xi = xi*(1 + xi_perc*(-1 + 2*rand));
    xi_ind = floor(N*xi);
    Je = sparse(N-1,1);
    Je(xi_ind) = k^2/(rho*A*h);
    
    % calculate output position
    xo_ind = floor(N*xo);
    c = sparse(1,N-1);
    c(xo_ind) = 1;
    
    % define fretboard
    Nc = N-1+Nfret;
    xfret = (L - L./(2.^((1:Nfret)/12))).';
    b = [Hback*ones(N-1,1); (Hfret+Hback)*ones(Nfret,1)];
    Gfret = sparse(N-1,Nfret);
    for m = 1:Nfret
        idx = floor(xfret(m)/h);
        w = xfret(m)/h - idx;
        Gfret(idx:idx+1,m) = [1-w;w];
    end
    Gc = (1/h)*[speye(N-1), Gfret];

    % initialise scheme variables
    u2 = zeros(N-1,1); % state
    u1 = u2;           % state
    u = u2;            % state
    y = zeros(Nf,1);   % output
    
    % calculate finger horizontal position
    gf = zeros(N-1,Nf);
    xfret = [0; xfret];
    hf = xfret(nf+1) - xfret(nf);
    xf = (xfret(nf) + 0.6*hf + xf_perc*hf*(-1+2*rand))/L;
    xf_ind = round(xf*N);
    gf(xf_ind,:) = 1/h;

    % define finger initial conditions
    uf2 = uf0;         % initial displacement
    uf1 = uf0 + k*vf0; % second displacement
    uf = uf1;
    
    % define finger force
    fin_st = fin_st*(1 + fin_st_perc*(-1 + 2*rand));
    f0amp = f0amp*(1 + f0amp_perc*(-1 + 2*rand));
    n_fin_st = floor(fin_st/k);
    n_fin_dur = floor(fin_dur/k);
    if n_fin_st + n_fin_dur > Nf
        n_fin_dur = Nf - n_fin_st;
    end
    f0 = [zeros(n_fin_st,1); f0amp*ones(n_fin_dur,1); zeros(Nf-n_fin_st-n_fin_dur,1)];
    Z = sparse(Nc+1, Nc+1);
    Z(end) = k^2/Mf;
    
    % define figure
    if plot_on == 1
        xax = (1:N-1)'*h; % x-axis for plotting
        fig_u = figure;
        hold on
        plot(xax, Hback*ones(size(xax)), 'b');
        for m = 2:Nfret+1
            plot([xfret(m), xfret(m)], [Hback, Hback+Hfret], 'b');
        end
        xlabel('x [m]', 'Interpreter', 'latex')
        ylabel('u [m]', 'Interpreter', 'latex')
        title("String state with fretboard (blue) and finger (green) collisions", 'Interpreter', 'latex')
        axis([0 L -0.005 0.005])
    end
    
    % main loop
    if vf0 == 0 && uf0 > 0
        n_st = min([min(find(f0 ~= 0)), min(find(fe ~= 0))]); % find excitation start
    else
        n_st = 1;
    end
    for n = n_st:Nf
        % calculate matrices for linear equation
        a = 0.5*k*sqrt(E*h/(rho*L))*Dxx*u1;
        A_inv = (1/(1+sig0*k))*(speye(N-1) - a*(a.')/(1+sig0*k+(a.')*a));
        C = (sig0*k-1)*speye(N-1) - a*(a.') - 2*sig1*k*Dxx;
        q = A_inv*(B*u1 + C*u2 + Je*fe(n));
        
        % collisions with fretboard and finger
        if calc_collisions == 1
			% finger horizontal position
			if n <= Nf-2
				mu_gf = 0.5*(gf(:,n) + gf(:,n+2));
			end
			if n == n_st
				mu_gf2 = mu_gf;
				mu_gf1 = mu_gf;
			end
			
			%
			G2 = [Gc, -mu_gf2];
			G1 = [Gc, -mu_gf1];
			G = [Gc, -mu_gf];
			
			%
			J = (k^2/(rho*A))*G1;
			J_ = A_inv*J;
			
			% calculate collision distances
			eta2 = [b; -uf2] - h*(G2.')*u2;
			eta_c2 = eta2(1:end-1);
			eta_f2 = eta2(end);
			eta1 = [b; -uf1] - h*(G1.')*u1;
			
			%
			P = sparse(Nc+1, Nc+1);
			P(Nc+1,Nc+1) = theta(eta1(end), Kf, af, bf)/(2*k);
			
			% calculate matrices for nonlinear equation
			M = Z + h*(G.')*J_;
			Q = speye(Nc+1) + M*P;
			I = -[zeros(Nc,1); -2*(uf1-uf2)+(k^2/Mf)*f0(n)] + h*(G.')*q - h*(G2.')*u2;
			
			% apply Newton-Raphson
			r = ones(Nc+1,1).*1e-3;
			iter = 0;
			step = ones(size(r));
			while (norm(step) > tol_nr) && (iter < max_iter)
				% implicit function matrices
				rc = r(1:end-1);
				rf = r(end);
				lambda = [(phi(eta_c2+rc,Kc,ac)-phi(eta_c2,Kc,ac))./rc;... 
				          (phi(eta_f2+rf,Kf,af)-phi(eta_f2,Kf,af))./rf];
				lambda(isnan(lambda) == 1) = 0;
				dlambda = [(dphi(eta_c2+rc,Kc,ac).*rc - (phi(eta_c2+rc,Kc,ac)-phi(eta_c2,Kc,ac)))./(rc.^2);... 
				           (dphi(eta_f2+rf,Kf,af).*rf - (phi(eta_f2+rf,Kf,af)-phi(eta_f2,Kf,af)))./(rf.^2)];
				dlambda(isnan(dlambda) == 1) = 0;
				F = Q*r + M*lambda + I;
				DF = Q + M*spdiags(dlambda(:),0,Nc+1,Nc+1);
					
				% Newton-Raphson step
				step = DF\F;
				r = r - step;
				iter = iter + 1;
			end
				
			% save result
			rc = r(1:end-1);
			rf = r(end);
			lambda = [(phi(eta_c2+rc,Kc,ac)-phi(eta_c2,Kc,ac))./rc;... 
			          (phi(eta_f2+rf,Kf,af)-phi(eta_f2,Kf,af))./rf];
			lambda(isnan(lambda) == 1) = 0;
			f = lambda + P*r;
			
			% update state, and insert current value of f
			u = q + J_*f;
			
			% update finger position
			uf = (k^2/Mf)*(f(end)-f0(n)) + 2*uf1 - uf2;
			
			% reset finger to initial position
			% (so it doesn't drift away without finger force)
			if uf > uf0
				uf = uf0;
				uf1 = uf;
				uf2 = uf;
			end
			
			% shift finger state
			mu_gf2 = mu_gf1;
			mu_gf1 = mu_gf;
			uf2 = uf1;
			uf1 = uf;
        else
            u = q;
        end
        
        % read output
        if otype == 1
            y(n) = c*u;
        elseif otype == 2
            y(n) = (c*u - c*u1)/k;
        end
    
        % plot
        if plot_on == 1
            % draw current state
            figure(fig_u)
            if exist('p','var')
                delete(p);
                delete(pf);
            end
            p = plot(xax, u, 'k');
            pf = plot(xax(xf_ind), uf, 'go');
            drawnow
        end
    
        % shift state
        u2 = u1;
        u1 = u;
    end
end

% knee function [1]
function y = kp(x)
    y = 0.5*(x + abs(x));
end

% collision potential function [1]
function y = phi(x, K, a)
    y = (K/(a+1))*kp(x).^(a+1);
end

% collision potential derivative [1]
function y = dphi(x, K, a)
    y = K*kp(x).^a;
end

% collision loss function [1]
function y = theta(x, K, a, b)
    y = K*b*kp(x)^a;
end

% dispersion relation approximate solution for wavenumber [2]
function out = xi_f(f, c, K)
    w = 2*pi*f;
    out = (-c^2 + sqrt(c^4 + 4*K^2*w^2)) / (2*K^2);
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