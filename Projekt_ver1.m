clear; clc; close all;

%% Parametry impulsu i próbkowania
tau = 1e-6;              % stała czasowa [s]
Fs = 20e6;               % częstość próbkowania [Hz]
Ts = 1/Fs;
tMax = 10*tau;
t = 0:Ts:tMax;
Nimp = length(t);

%% Generacja impulsu H(t) (kształtowanie)
H_t = (t./tau.^2).*exp(-t/tau);   % odpowiedź systemu
H_t = H_t / max(H_t);             % normalizacja do max = 1 (czyli szukana amplituda)

%% Próbkowanie impulsu z losowym przesunięciem fazy
delay = rand()*Ts;                % losowe przesunięcie
t_sampled = delay:Ts:(tMax+delay);
H_sampled = (t_sampled./tau.^2).*exp(-t_sampled/tau);
H_sampled = H_sampled / max(H_sampled);

% Zakładamy że mamy tylko kilka próbek z impulsu
sampledSegment = H_sampled(1:10);

%% Projektowanie filtru FIR odwrotnego do H(f)
Nf = 2*31+1;
K = 2*7+1;
window = hamming(Nf)';

% Tworzymy H(f) jako prostokąt (to uproszczenie dla przykładu)
zerosBefore = floor((Nf-K)/2);
zerosAfter  = ceil((Nf-K)/2);
requestedH = [zeros(1,zerosBefore), ones(1,K), zeros(1,zerosAfter)];

% Symetryczne indeksy
indexesSym = -ceil(Nf/2)+1 : floor(Nf/2);

figure;
subplot(5,1,1); stem(indexesSym, requestedH); title('Zadana funkcja H(f)'); xlim([indexesSym(1), indexesSym(end)]);

% Shift do zakresu (0, fs)
shiftedH = ifftshift(requestedH);
subplot(5,1,2); stem(0:Nf-1, shiftedH); title('Przesunięta H(f)');

% IFFT -> impuls odwrotny
timeResponse = ifft(shiftedH);
FIRcoefs = fftshift(timeResponse); % współczynniki filtru
FIRcoefs = FIRcoefs .* window;     % okno Hamminga
subplot(5,1,3); stem(indexesSym, timeResponse); title('Odwrotna odpowiedź impulsowa (IFFT)');
subplot(5,1,4); stem(0:Nf-1, FIRcoefs); title('Filtr FIR odwrotny');

% Sprawdzenie przenoszenia
realH = fft(FIRcoefs, 1000);
subplot(5,1,5); plot(20*log10(abs(realH))); title('Charakterystyka filtru [dB]');

%% Filtracja próbek (rekonstrukcja amplitudy)
% Dopełniamy próbkowany impuls zerami do długości filtru
x = [sampledSegment, zeros(1, Nf-length(sampledSegment))];
y = conv(x, FIRcoefs, 'same'); % filtracja
amp_rec = max(y);              % rekonstrukcja amplitudy

%% Porównanie
fprintf('Prawdziwa amplituda: %.4f\n', 1.0);
fprintf('Zrekonstruowana amplituda: %.4f\n', amp_rec);

% Wykres porównawczy
figure;
plot(x, '-o'); hold on;
plot(y, '-x');
legend('Oryginalne próbki','Po filtracji');
title('Rekonstrukcja amplitudy z próbek');
xlabel('Próbka'); ylabel('Amplituda');
