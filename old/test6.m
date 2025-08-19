clear; clc; close all;

%% 1. Parametry
fs = 1000;             % Częstotliwość próbkowania [Hz]
T = 1/fs;
t = 0:T:0.01;          % 10 ms
tau = 0.001;           % Stała czasowa shaping'u
n = 2;                 % Parametr kształtu impulsu
true_amplitude = 1;    % Szukana amplituda

%% 2. Odpowiedź analogowa H(t)
Ht = (t.^n).*exp(-t/tau);    % Quasi-gaussowski impuls
Ht = Ht / max(Ht) * true_amplitude;  % Normalizacja do amplitudy
L = length(Ht);

figure;
subplot(4,1,1);
plot(t*1000, Ht);
title('Odpowiedź systemu H(t)');
xlabel('t [ms]');
ylabel('Amplituda');

%% 3. Próbkowanie w losowej fazie
sampling_interval = 5;  % co ile próbek próbkujemy (niższy fs ADC)
offset = randi([1, sampling_interval]);  % losowa faza
sampled = Ht(offset:sampling_interval:end);
t_sampled = t(offset:sampling_interval:end);

subplot(4,1,2);
stem(t_sampled*1000, sampled, 'filled');
title('Próbkowany sygnał (losowa faza)');
xlabel('t [ms]');
ylabel('Amplituda');

%% 4. Projektowanie filtru odwrotnego FIR
N = 2*31+1;   % Rząd filtru (nieparzysty)
Hf = fft(Ht, N);         % Odpowiedź częstotliwościowa
Hf = Hf + 1e-6;          % Uniknięcie dzielenia przez 0
invHf = 1 ./ Hf;         % Odwrotna funkcja przenoszenia

% Transformata odwrotna
invHt = real(ifft(invHf));
invHt = fftshift(invHt);        % Przesunięcie zer na środek

% Okno Hamming'a
window = hamming(N)';
invHt_windowed = invHt .* window;

subplot(4,1,3);
stem(0:N-1, invHt_windowed, 'filled');
title('Współczynniki filtru FIR (po oknie)');
xlabel('Indeks próbki');

%% 5. Filtracja (rekonstrukcja amplitudy)
reconstructed = conv(sampled, invHt_windowed, 'same');

subplot(4,1,4);
plot(t_sampled*1000, reconstructed, '-o');
title('Zrekonstruowany sygnał');
xlabel('t [ms]');
ylabel('Amplituda');

%% 6. Ocena rekonstrukcji
[recovered_amplitude, idx] = max(reconstructed);
fprintf('---\n');
fprintf('Prawdziwa amplituda:        %.4f\n', true_amplitude);
fprintf('Zrekonstruowana amplituda:  %.4f\n', recovered_amplitude);
fprintf('Błąd względny:              %.2f %%\n', abs((recovered_amplitude - true_amplitude)/true_amplitude)*100);
fprintf('---\n');
