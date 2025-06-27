% Parametry modelu
A = 1;         % amplituda
tau = 1e-6;    % stała czasowa [s]
n = 1;         % rząd shapera

% Oś czasu
t_max = 10 * tau;
dt = 1e-8;
t = 0:dt:t_max;

% H(t)
H = A * (t / tau).^n .* exp(-t / tau);
H(t < 0) = 0;

% Normalizacja amplitudy
H = H / max(H);

figure;
plot(t, H);
xlabel('Czas [s]'); ylabel('H(t)'); title('Odpowiedź impulsowa shapera');
% Modelujemy impuls delta w chwili t = 0
d = zeros(size(t));
d(1) = 1/dt;  % aproksymacja delta funkcji

% Splot
Vobs = conv(d, H, 'same');

% Dodajmy szum (opcjonalnie)
noise_amp = 0.01;
Vobs_noisy = Vobs + noise_amp * randn(size(Vobs));

figure;
plot(t, Vobs_noisy);
xlabel('Czas [s]'); ylabel('V_{obs}(t)'); title('Sygnał obserwowany z szumem');
fs = 10e6;                   % częstotliwość próbkowania [Hz]
Ts = 1/fs;
t0 = rand * Ts;              % losowy offset fazowy
ts = t0:Ts:t_max;            % chwile próbkowania

% Interpolacja próbek
V_samples = interp1(t, Vobs_noisy, ts, 'linear', 0);

figure;
plot(t, Vobs_noisy); hold on;
stem(ts, V_samples, 'r');
xlabel('Czas [s]'); ylabel('Sygnał');
title('Próbkowanie z losowym offsetem fazowym');
legend('V_{obs}(t)', 'V[n]');
N = 64;  % rząd filtru FIR

% FFT H(t)
H_fft = fft(H, 2^nextpow2(length(H)));
f = linspace(0, fs, length(H_fft));

% Idealna odpowiedź filtru dekonwolucyjnego
H_fft_inv = 1 ./ H_fft;
H_fft_inv(abs(H_fft) < 1e-3) = 0;  % maskowanie niestabilnych częstotliwości

% Filtr FIR przez IFFT (i okno Hamminga)
h_inv = real(ifft(H_fft_inv));
h_inv = h_inv(1:N) .* hamming(N)';

% Zastosowanie filtru FIR
V_rec = conv(V_samples, h_inv, 'same');

figure;
stem(V_rec);
xlabel('n'); ylabel('d_{rec}[n]');
title('Zrekonstruowany sygnał po dekonwolucji');
amp_est = max(V_rec);
fprintf('Oszacowana amplituda: %.3f\n', amp_est);
