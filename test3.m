clear; close all; clc;

%% Parametry
T = 25e-9;          % okres próbkowania 25 ns (40 MHz)
t_end = 1e-6;       % długość obserwacji 1 us
t0 = 0;
t = linspace(t0, t_end, 2^14); % gęsta siatka czasowa do "ciągłego" sygnału

% Odpowiedź impulsowa detektora (np. funkcja eksponencjalna + oscylacje)
tau = 100e-9; % stała czasowa
f_osc = 5e6;  % 5 MHz oscylacje
H = (t >= 0) .* exp(-t/tau) .* sin(2*pi*f_osc*t);

% FFT i projektowanie filtru odwrotnego FIR
Nfir = 101;
Nfft = 2^nextpow2(length(H)*2);
Hf = fft(H, Nfft);

alpha = 1e-3; % regularyzacja
invHf = conj(Hf) ./ (abs(Hf).^2 + alpha);
h_inv_full = real(ifft(invHf));
h_inv = fftshift(h_inv_full);

mid = floor(length(h_inv)/2)+1;
fir_coeffs = h_inv(mid - floor(Nfir/2) : mid + floor(Nfir/2));

% Pokazujemy współczynniki filtru
figure;
plot(fir_coeffs, '-o');
title('Współczynniki filtru FIR');
grid on;

%% Symulacja amplitud i rekonstrukcji
num_trials = 1000;
Qdet_vals = 5 + rand(1, num_trials); % wartości Qdet losowe z przedziału [5,6]
rec_amplitudes = zeros(1, num_trials);

noise_level = 0.01;

for i = 1:num_trials
    Qdet = Qdet_vals(i);

    % Sygnal ciągły
    sig = Qdet * H;

    % Dodanie szumu
    sig_noisy = sig + noise_level * randn(size(sig));

    % Losowa faza próbkowania
    phase_shift = rand()*T;
    t_sampled = t0 + phase_shift : T : t_end;

    % Próbkowanie
    sampled_sig = interp1(t, sig_noisy, t_sampled, 'linear', 0);

    % Rekonstrukcja przez filtrację
    recon_sig_full = conv(sampled_sig, fir_coeffs, 'full');
    recon_sig = recon_sig_full(floor(Nfir/2)+1 : floor(Nfir/2)+length(sampled_sig));

    % Maksimum zrekonstruowanego sygnału = estymowana amplituda
    rec_amplitudes(i) = max(recon_sig);
end

%% Wykres Qdet vs rekonstruowane amplitudy
figure;
plot(Qdet_vals, rec_amplitudes, '.');
xlabel('Q_{det} (prawdziwa amplituda)');
ylabel('Zrekonstruowana amplituda');
title('Rekonstrukcja amplitudy z sygnału próbkowanego');
grid on;

% Dopasowanie prostej regresji
p = polyfit(Qdet_vals, rec_amplitudes, 1);
hold on;
plot(Qdet_vals, polyval(p, Qdet_vals), 'r', 'LineWidth', 2);
legend('Dane', sprintf('Regresja: y = %.2fx + %.2f', p(1), p(2)));

figure;
plot(t_sampled*1e6, recon_sig);
title('Sygnał po rekonstrukcji');
xlabel('Czas [\mus]');
ylabel('Amplituda');
grid on;

figure;
plot(t*1e6, sig_noisy, 'b', 'DisplayName', 'Sygnał z szumem');
hold on;
stem(t_sampled*1e6, sampled_sig, 'r', 'DisplayName', 'Próbkowanie');
legend;
xlabel('Czas [\mus]');
ylabel('Amplituda');
title('Sygnał i próbki');
grid on;

figure;
subplot(3,1,1);
plot(t*1e6, H);
title('Odpowiedź impulsowa H(t) - shaper');
xlabel('Czas [\mus]');
ylabel('Amplituda [arb. units]');
grid on;

% subplot(3,1,2);
% plot(t*1e6, signal);
% title('Modelowany sygnał wyjściowy (splot \delta z H)');
% xlabel('Czas [\mus]');
% ylabel('Amplituda');
% grid on;
% 
% subplot(3,1,3);
% plot(Qdet_values, reconstructed_amplitudes, 'o-', 'LineWidth', 2);
% hold on;
% plot(Qdet_values, Qdet_values, 'k--', 'LineWidth', 1.5);
% title('Rekonstrukcja amplitudy sygnału');
% xlabel('Amplituda wejściowa Q_{det}');
% ylabel('Odzyskana amplituda');
% legend('Odzyskana', 'Idealna', 'Location', 'NorthWest');
% grid on;
% 
% figure;
% plot(Qdet_values, errors, 'r*-');
% title('Błąd rekonstrukcji amplitudy');
% xlabel('Amplituda wejściowa Q_{det}');
% ylabel('Błąd bezwzględny');
% grid on;