% Parametry sygnału
A = 1;          % Skalujący czynnik
tau = 0.1;      % Stała czasowa
n = 1;          % Rząd shapera (np. n=1 dla CR-RC)

% Częstotliwość próbkowania (większa niż fNyquist)
fs = 2 / tau * 5;   % Przykładowa częstotliwość próbkowania
Ts = 1/fs;

% Dostosowanie zakresu czasowego
t_start = -0.1;
t_end = 0.5;
t = t_start : Ts : t_end;

% Definicja funkcji H(t)
H_t = A * (t ./ tau).^n .* exp(-t / tau);

% Generacja sygnału d(t) – np. impuls jednostkowy
d_t = zeros(size(t));
d_t(ceil(length(t)/2)) = 1; % Impuls w środku

% Obliczanie splotu Vobs = d(t) * H(t)
V_obs_t = conv(d_t, H_t);

% Przygotowanie do próbkowania – wybranie odpowiednich próbek
t_idx = 1 : length(V_obs_t);
t_sampled = t(1:end); % Próby w czasie

% Przekształcenie Fouriera sygnałów (np. dla dekonwolucji)
f_fft = (-length(V_obs_t)/2 : length(V_obs_t)/2 - 1) * fs / length(V_obs_t);
V_fft = fft(V_obs_t);
H_fft = fft(H_t);

% Zastosowanie filtru dekonwolucyjnego: D(f) = Vobs(f) / H(f)
D_rec_fft = V_fft ./ H_fft;

% Rekonstrukcja sygnału d(t)
d_rec_t = ifft(D_rec_fft);

% Wyświetlenie wyników
figure;
subplot(3, 1, 1);
plot(t_sampled, real(V_obs_t));
title('Obserwowany sygnał V_obs(t)');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(3, 1, 2);
plot(f_fft, abs(H_fft), 'r', 'LineWidth', 2);
title('Widmo H(f)');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda');

subplot(3, 1, 3);
plot(t_sampled, real(d_rec_t));
title('Rekonstruowany sygnał d_rec(t)');
xlabel('Czas [s]');
ylabel('Amplituda');
