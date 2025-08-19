clear; clc; close all;

%% Parametry impulsu i próbkowania
A = 1;              % Współczynnik skalujący
n = 2;              % Wykładnik potęgi
tau = 1e-6;         % Stała czasowa [s]
Fs = 20e6;          % Częstość próbkowania [Hz]
Ts = 1/Fs;
tMax = 10*tau;
t = 0:Ts:tMax;

%% Generacja impulsu H(t)
H_t = A * (t ./ tau).^n .* exp(-t / tau);
H_t = H_t / max(H_t);  % Normalizacja (maksimum = 1)

% Wykres funkcji H(t)
figure;
plot(t, H_t);
xlabel('Czas t [s]');
ylabel('H(t)');
title('Funkcja H(t) = A*(t/tau)^n * exp(-t/tau)');
grid on;

%% Próbkowanie impulsu z losowym przesunięciem fazy
delay = rand() * Ts;                         % Losowy offset fazy
t_sampled = delay : Ts : (tMax + delay);     % Przesunięta siatka próbkowania

% Próbkowanie oryginalnego impulsu H(t)
H_sampled = A * (t_sampled ./ tau).^n .* exp(-t_sampled / tau);
H_sampled = H_sampled / max(H_sampled);      % Normalizacja

% Zakładamy tylko kilka pierwszych próbek z impulsu
num_samples = 10;
sampledSegment = H_sampled(1:num_samples);

%% Projektowanie filtru FIR odwrotnego do H(f)
Nf = 2*31+1;         % Długość filtru (nieparzysta)
K = 2*7+1;           % Szerokość pasma (ilość jedynek w "prostokącie")
window = hamming(Nf)'; % Okno Hamminga

% Zadana funkcja H(f) (prostokąt)
zerosBefore = floor((Nf - K)/2);
zerosAfter  = ceil((Nf - K)/2);
requestedH = [zeros(1,zerosBefore), ones(1,K), zeros(1,zerosAfter)];

% Symetryczne indeksy
indexesSym = -ceil(Nf/2)+1 : floor(Nf/2);

% Wykresy pośrednie
figure;
subplot(5,1,1);
stem(indexesSym, requestedH, 'filled');
title('Zadana funkcja H(f)'); xlim([indexesSym(1), indexesSym(end)]);

% Shift do zakresu (0, fs)
shiftedH = ifftshift(requestedH);
subplot(5,1,2);
stem(0:Nf-1, shiftedH, 'filled');
title('Przesunięta H(f)');

% IFFT -> impuls odwrotny
timeResponse = ifft(shiftedH, 'symmetric');  % wymuszamy rzeczywiste wartości
FIRcoefs = fftshift(timeResponse);           % Środek na środku
FIRcoefs = FIRcoefs .* window;               % Zastosowanie okna
subplot(5,1,3);
stem(indexesSym, timeResponse, 'filled');
title('Odwrotna odpowiedź impulsowa (IFFT)');

subplot(5,1,4);
stem(0:Nf-1, FIRcoefs, 'filled');
title('Filtr FIR odwrotny');

% Sprawdzenie charakterystyki częstotliwościowej filtru
realH = fft(FIRcoefs, 1000); % interpolowana DFT
subplot(5,1,5);
plot(linspace(0, Fs/2, 500), 20*log10(abs(realH(1:500))));
title('Charakterystyka filtru [dB]');
xlabel('Częstotliwość [Hz]');
ylabel('Wzmocnienie [dB]'); grid on;

%% Filtracja próbek (rekonstrukcja amplitudy)
% Dopełnienie zerami, jeśli potrzeba
x = [sampledSegment, zeros(1, Nf - length(sampledSegment))];

% Filtracja
y = conv(x, FIRcoefs, 'same');
amp_rec = max(y);  % Rekonstrukcja amplitudy

%% Porównanie amplitud
fprintf('Prawdziwa amplituda:       %.4f\n', 1.0);
fprintf('Zrekonstruowana amplituda: %.4f\n', amp_rec);
fprintf('Błąd względny:             %.2f%%\n', 100 * abs(amp_rec - 1.0));

% Wykres porównawczy
figure;
plot(x, '-o'); hold on;
plot(y, '-x');
legend('Oryginalne próbki','Po filtracji');
title('Rekonstrukcja amplitudy z próbek');
xlabel('Numer próbki');
ylabel('Amplituda');
grid on;
