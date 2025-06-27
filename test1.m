clear; clc; close all;

%% PARAMETRY SYSTEMU I SYMULACJI
A = 1;                  % amplituda
n = 2;                  % rząd shaping function
tau = 1e-6;             % stała czasowa
Tsim = 10 * tau;        % długość symulacji

% Analogowe próbkowanie sygnału (wysoka rozdzielczość)
Fs_analog = 1e9;
dt = 1 / Fs_analog;
t = 0:dt:Tsim;

% Transfer function H(t)
H_t = A * (t ./ tau).^n .* exp(-t / tau);
H_t(t < 0) = 0;
H_t = H_t / max(H_t);   % normalizacja do 1

%% TRANSFORMATA FOURIERA H(f)
H_f = fft(H_t);
f = (0:length(H_t)-1) * Fs_analog / length(H_t);

%% FILTR WIENERA (DEKONWOLUCJA)
% Parametry szumu
noise_amplitude_rms = 0.01;
k = noise_amplitude_rms^2;

% Odwrotność H(f) z regularyzacją
H_f_inv = conj(H_f) ./ (abs(H_f).^2 + k);

% Impuls filtru (długa odpowiedź)
h_inv_long = real(ifft(H_f_inv));
h_inv_long = fftshift(h_inv_long);  % centrowanie

% Przycinanie i okno Hamming'a
Nf = 63;  % długość filtra FIR (nieparzysta)
center = floor(length(h_inv_long)/2);
h_short = h_inv_long(center - floor(Nf/2) : center + floor(Nf/2));
FIRcoefs = h_short .* hamming(Nf)';

%% PRÓBKOWANIE I ODTWORZENIE AMPLITUDY
Fs = 20e6;              % częstotliwość ADC
Ts = 1 / Fs;
delay = rand() * Ts;    % losowy offset fazy

% Punkty próbkowania
t_sampled = delay : Ts : (t(end) + delay);
H_sampled = A * (t_sampled ./ tau).^n .* exp(-t_sampled / tau);
H_sampled(t_sampled < 0) = 0;
H_sampled = H_sampled / max(H_sampled);  % normalizacja

% Dodaj szum i weź fragment
num_samples = 10;
x = H_sampled(1:num_samples) + noise_amplitude_rms * randn(1, num_samples);

% Zero-padding i filtracja
x_padded = [x, zeros(1, Nf - length(x))];
y = conv(x_padded, FIRcoefs, 'same');

% Rekonstrukcja amplitudy jako maksimum
amp_rec = max(y);

%% WYNIKI
fprintf('--- WYNIKI REKONSTRUKCJI ---\n');
fprintf('Prawdziwa amplituda:       %.4f\n', A);
fprintf('Zrekonstruowana amplituda: %.4f\n', amp_rec);
fprintf('Błąd względny:             %.2f%%\n', 100 * abs(amp_rec - A) / A);

%% WIZUALIZACJA
figure;
subplot(2,1,1);
plot(t*1e6, H_t); grid on;
xlabel('Czas [\mus]'); ylabel('H(t)');
title('Shaping function H(t)');

subplot(2,1,2);
plot(f(1:floor(end/2))/1e6, abs(H_f(1:floor(end/2))));
xlabel('Częstotliwość [MHz]'); ylabel('|H(f)|');
title('Moduł transformaty Fouriera H(f)');
grid on;

figure;
stem(FIRcoefs, 'filled');
xlabel('Numer próbki'); ylabel('Amplituda');
title('Filtr FIR - odpowiedź impulsowa (Wiener)');
grid on;

figure;
plot(x, '-o'); hold on;
plot(y, '-x');
legend('Próbki zaszumione', 'Po filtracji');
xlabel('Numer próbki'); ylabel('Amplituda');
title('Rekonstrukcja amplitudy');
grid on;
