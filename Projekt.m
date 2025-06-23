%% 1. Definicja parametrów symulacji

clear; close all; clc;

% Czas symulacji
T_sim = 200e-9; % Czas trwania symulacji w sekundach (np. 200 ns)
Fs_analog = 10e9; % Częstotliwość próbkowania dla "ciągłego" sygnału (bardzo wysoka, np. 10 GHz)
dt_analog = 1/Fs_analog;
t_analog = 0:dt_analog:T_sim-dt_analog; % Wektor czasu dla sygnału "ciągłego"

% Parametry sygnału z detektora (delta Kroneckera)
amplitudes = [0.1, 0.5, 1, 2, 5, 10]; % Przykładowe amplitudy ładunku (np. w pC)
num_events_per_amplitude = 100; % Liczba zdarzeń do symulacji dla każdej amplitudy

% Parametry ADC
Fs_adc_range = [50e6, 100e6, 200e6, 500e6]; % Zakres częstotliwości próbkowania ADC (np. od 50 MHz do 500 MHz)
% Dla uproszczenia na początku skupimy się na jednej częstotliwości ADC,
% a potem rozszerzymy analizę na zakres.
Fs_adc_selected = Fs_adc_range(2); % Wybierz jedną do początkowej symulacji

% Parametry szumu
noise_amplitude_rms = 0.01; % RMS szumu (np. w Voltach)


%% 2. Modelowanie Funkcji Przenoszenia Elektroniki Odczytu H(t)
% To jest KLUCZOWY element, który musisz dokładnie określić.
% Poniżej przykład prostej funkcji quasi-gaussowskiej (np. z kaskady CR-RC^n)

% Funkcja CR-RC^n jest często używana do modelowania.
% Przykład dla CR-RC (n=1) lub CR-RC^2 (n=2)
% H(t) = C * (t/tau)^n * exp(-t/tau) dla t >= 0
% Gdzie C jest stałą normalizacyjną, tau to stała czasowa shaping-u.

tau = 20e-9; % Stała czasowa shapera (np. 20 ns)
n_shaper = 2; % Rząd shapera (np. 2 dla CR-RC^2, co daje kształt bardziej zbliżony do Gaussa)

% Definicja H(t)
H_t = zeros(size(t_analog));
valid_indices = t_analog >= 0; % Impuls zaczyna się w t=0
H_t(valid_indices) = (t_analog(valid_indices) / tau).^n_shaper .* exp(-t_analog(valid_indices) / tau);

% Normalizacja H(t) tak, aby maksimum było równe 1 dla danego ładunku (po wzmocnieniu)
% W rzeczywistości amplituda zależy od wzmocnienia, ale dla modelowania możemy ją unormować
% do umownej "jednostki" i potem skalować przez "amplitudy"
[~, max_idx] = max(H_t);
H_t = H_t / H_t(max_idx); % Normalizacja do 1

% Wizualizacja H(t)
figure;
plot(t_analog * 1e9, H_t);
title('Funkcja przenoszenia elektroniki odczytu H(t)');
xlabel('Czas [ns]');
ylabel('Amplituda (unormowana)');
grid on;

% Obliczenie Transformacji Fouriera H(f) dla celów projektowania filtra
% Używamy FFT dla dyskretnego sygnału
H_f_analog = fft(H_t);
f_analog = (0:length(H_t)-1) * Fs_analog / length(H_t); % Wektor częstotliwości

% Wizualizacja |H(f)|
figure;
plot(f_analog(1:floor(length(f_analog)/2))/1e6, abs(H_f_analog(1:floor(length(H_f_analog)/2))));
title('Moduł transformaty Fouriera H(f)');
xlabel('Częstotliwość [MHz]');
ylabel('|H(f)|');
grid on;


%% 6. Projektowanie Filtra Dekonwolucyjnego 1/H(f) (FIR)

% Długość filtra FIR - ma wpływ na dokładność i złożoność
% Dłuższy filtr = lepsze odwzorowanie 1/H(f), ale większe opóźnienie i złożoność
fir_len = 50; % Przykładowa długość filtra FIR

% !!!!!!! KLUCZOWY PUNKT: Projektowanie filtra 1/H(f) !!!!!!!
% Bezpośrednie 1./H_f_analog prowadzi do problemów z szumem i niestabilności.
% Musisz zastosować regularyzację lub bardziej zaawansowane techniki.
% Poniżej prosty przykład z progowaniem (truncation) lub filtrem Wienera.

% Metoda 1: Proste odwrócenie z progowaniem małych wartości (NIEZALECANE, ale dla ilustracji)
% H_f_inv = zeros(size(H_f_analog));
% threshold = max(abs(H_f_analog)) * 0.01; % Próg dla małych wartości
% valid_freq_idx = abs(H_f_analog) > threshold;
% H_f_inv(valid_freq_idx) = 1 ./ H_f_analog(valid_freq_idx);

% Metoda 2: Filtr Wienera (bardzo często używana w dekonwolucji)
% G(f) = H*(f) / (|H(f)|^2 + S_n(f)/S_x(f))
% Gdzie S_n(f) to widmowa gęstość mocy szumu, S_x(f) to widmowa gęstość mocy sygnału
% Dla uproszczenia, możemy użyć stałej k jako stosunku szumu do sygnału:
% G(f) = H*(f) / (|H(f)|^2 + k)
% Musisz określić k (power_spectral_density_ratio) na podstawie szumu i sygnału.
% Możemy przyjąć uproszczoną wersję gdzie k jest stałą regularyzacji.

regularization_k = (noise_amplitude_rms)^2 / (max(amplitudes))^2; % Przykład, dostosuj!

H_f_inv_wiener = conj(H_f_analog) ./ (abs(H_f_analog).^2 + regularization_k);

% Przejście do dziedziny czasu - odpowiedź impulsowa filtra FIR
% Używamy IFFT, a następnie skracamy do długości fir_len.
% W tym miejscu możesz również użyć funkcji do projektowania filtrów FIR, np. `fir1`
% ale to wymagałoby określenia pasm przenoszenia/tłumienia.
% Bezpośrednie przejście z dziedziny częstotliwości:

impulse_response_fir_long = ifft(H_f_inv_wiener);
% Centrowanie i skrócenie odpowiedzi impulsowej
% Zwykle chcemy, aby impuls był symetryczny wokół środka dla liniowej fazy
impulse_response_fir = fftshift(impulse_response_fir_long); % Przesunięcie zera częstotliwości
impulse_response_fir = impulse_response_fir(floor(end/2) - floor(fir_len/2) + 1 : floor(end/2) + ceil(fir_len/2));

% Normalizacja filtru (opcjonalnie, ale może pomóc w zachowaniu amplitudy)
impulse_response_fir = real(impulse_response_fir); % Pozbywamy się małych części urojonych

% Wizualizacja odpowiedzi impulsowej filtra FIR
figure;
stem(impulse_response_fir);
title('Odpowiedź impulsowa filtru FIR');
xlabel('Numer próbki');
ylabel('Wartość');
grid on;


%% 3-5, 7. Główna pętla symulacyjna: Generacja, próbkowanie, szum, rekonstrukcja
% Będziemy iterować po różnych amplitudach i częstotliwościach próbkowania ADC

reconstruction_errors_amplitude = zeros(length(amplitudes), 1);
reconstruction_errors_sampling = zeros(length(Fs_adc_range), 1);

% Dla uproszczenia, najpierw analizujemy tylko dla jednej Fs_adc_selected
fprintf('\n--- Symulacja dla Fs_ADC = %.1f MHz ---\n', Fs_adc_selected / 1e6);

mean_reconstructed_amplitudes = zeros(length(amplitudes), 1);
std_reconstructed_amplitudes = zeros(length(amplitudes), 1);

for amp_idx = 1:length(amplitudes)
    current_amplitude = amplitudes(amp_idx);
    reconstructed_amplitudes_for_current_amp = zeros(num_events_per_amplitude, 1);

    for event_idx = 1:num_events_per_amplitude
        % Generacja sygnału "ciągłego" (splot delty z H(t))
        % W idealnym przypadku to po prostu H_t skalowane przez amplitudę
        % Ale jeśli delta byłaby w innym miejscu, należałoby zrobić splot.
        % Tutaj zakładamy, że delta jest w t=0.
        continuous_signal = current_amplitude * H_t;

        % Próbkowanie ADC w losowej fazie
        % Losowe przesunięcie początkowe
        max_shift_samples_analog = floor(Fs_analog / Fs_adc_selected); % Max przesunięcie w próbkach "analogowych"
        random_shift_samples_analog = randi([0, max_shift_samples_analog - 1]);
        
        t_start_adc = t_analog(random_shift_samples_analog + 1);
        t_adc = t_start_adc : (1/Fs_adc_selected) : T_sim;
        
        % Interpolacja sygnału "ciągłego" na punkty próbkowania ADC
        % Użyjemy interp1 dla interpolacji liniowej, lub nearest dla najbliższego sąsiada
        sampled_signal_clean = interp1(t_analog, continuous_signal, t_adc, 'linear', 0);
        
        % Dodanie szumu Gaussa
        sampled_signal_noisy = sampled_signal_clean + noise_amplitude_rms * randn(size(sampled_signal_clean));

        % !!!!!!!! KLUCZOWY PUNKT: Zastosowanie filtru FIR do dekonwolucji !!!!!!!!
        % Funkcja `conv` wykonuje splot.
        % Wynik splotu będzie dłuższy niż sygnał wejściowy.
        % Musimy wybrać odpowiedni fragment, który odpowiada dekonwolucji.
        
        reconstructed_signal_raw = conv(sampled_signal_noisy, impulse_response_fir);
        
        % Znalezienie zrekonstruowanej amplitudy
        % Idealnie, po dekonwolucji powinniśmy otrzymać impuls Diraca (lub jego przybliżenie).
        % Amplituda tego impulsu będzie naszą zrekonstruowaną amplitudą.
        % Szukamy maksimum lub wartości w przewidywanym miejscu impulsu.
        
        % Problem: filtr FIR wprowadza opóźnienie. Należy to uwzględnić.
        % Opóźnienie = (length(impulse_response_fir) - 1) / 2 próbki (dla filtra liniowofazowego)
        
        delay_samples = floor(length(impulse_response_fir) / 2);
        
        % Spróbujmy znaleźć maksimum w pewnym oknie czasowym,
        % lub po prostu weźmy wartość w miejscu, gdzie spodziewamy się delty.
        
        % Znalezienie przybliżonej pozycji "zrekonstruowanej delty"
        % Oryginalny impuls był w okolicy t=0.
        % sampled_signal_noisy zaczyna się od t_start_adc.
        % Po konwolucji z filtrem, szczyt będzie przesunięty.
        
        % W prostym przypadku, jeśli filtrujemy cały sygnał i oczekujemy jeden impuls,
        % możemy szukać globalnego maksimum po pewnym opóźnieniu.
        
        % Przybliżona pozycja szczytu zrekonstruowanego impulsu
        % Jest to bardzo wrażliwy punkt i wymaga precyzyjnego ustawienia.
        % Można estymować szczyt zrekonstruowanego impulsu.
        
        % Dla uproszczenia: weźmy maksimum zrekonstruowanego sygnału po pewnym opóźnieniu.
        % To nie jest idealne, ale daje wstępne wyniki.
        
        [peak_val, peak_idx_raw] = max(reconstructed_signal_raw);
        % Należy uważać na artefakty na początku i końcu.
        % Prawdziwy szczyt powinien być w okolicy t_start_adc + delay_samples / Fs_adc_selected.
        
        % Możliwe jest również dopasowanie do znanej postaci impulsu Diraca po dekonwolucji.
        % Albo, jeśli filtr jest dobrze zaprojektowany, to szczyt 'reconstructed_signal_raw'
        % powinien być blisko faktycznej amplitudy.
        
        % Bardziej zaawansowana metoda: fitowanie impulsu lub szukanie maksimum w oknie
        
        % Na razie weźmy po prostu maksimum (po odrzuceniu początkowych/końcowych artefaktów)
        % Opóźnienie wprowadzone przez filtr FIR powoduje, że szczyt pojawi się później.
        % Sygnał po dekonwolucji powinien wyglądać jak krótki impuls.
        
        % Zakładamy, że szczyt jest w środku użytecznej części sygnału.
        % Długość sygnału po konwolucji: N_in + N_fir - 1
        % Wybieramy fragment odpowiadający oryginalnemu sygnałowi.
        
        reconstructed_signal_aligned = reconstructed_signal_raw(delay_samples + 1 : end - delay_samples);
        if isempty(reconstructed_signal_aligned)
            reconstructed_amplitudes_for_current_amp(event_idx) = 0; % W przypadku błędu
        else
            reconstructed_amplitudes_for_current_amp(event_idx) = max(reconstructed_signal_aligned);
        end
        
    end
    
    mean_reconstructed_amplitudes(amp_idx) = mean(reconstructed_amplitudes_for_current_amp);
    std_reconstructed_amplitudes(amp_idx) = std(reconstructed_amplitudes_for_current_amp);
    
    % Błąd rekonstrukcji dla danej amplitudy (np. odchylenie standardowe zrekonstruowanych wartości)
    % Lub średni błąd bezwzględny / kwadratowy
    reconstruction_errors_amplitude(amp_idx) = std(reconstructed_amplitudes_for_current_amp - current_amplitude);
    fprintf('  Amplituda: %.2f, Średnia zrekonstruowana: %.2f, Std dev: %.4f, Błąd: %.4f\n', ...
            current_amplitude, mean_reconstructed_amplitudes(amp_idx), ...
            std_reconstructed_amplitudes(amp_idx), reconstruction_errors_amplitude(amp_idx));
end

%% 8. Analiza Dokładności Rekonstrukcji w funkcji amplitudy
figure;
errorbar(amplitudes, mean_reconstructed_amplitudes, std_reconstructed_amplitudes, 'o-');
hold on;
plot(amplitudes, amplitudes, 'k--', 'DisplayName', 'Idealna rekonstrukcja');
title('Rekonstruowana Amplituda vs. Rzeczywista Amplituda');
xlabel('Rzeczywista Amplituda');
ylabel('Zrekonstruowana Amplituda (średnia +/- std dev)');
legend('show');
grid on;

figure;
plot(amplitudes, reconstruction_errors_amplitude, 'o-');
title('Błąd Rekonstrukcji (Std Dev) vs. Amplituda');
xlabel('Rzeczywista Amplituda');
ylabel('Błąd Rekonstrukcji (Std Dev)');
grid on;


% --- Rozszerzenie analizy na różne częstotliwości próbkowania ADC ---
fprintf('\n--- Analiza w funkcji częstotliwości próbkowania ADC ---\n');

% Wybieramy jedną, reprezentatywną amplitudę do tej analizy
representative_amplitude = amplitudes(end); % Np. największa amplituda

reconstruction_errors_sampling = zeros(length(Fs_adc_range), 1);
mean_reconstructed_sampling = zeros(length(Fs_adc_range), 1);
std_reconstructed_sampling = zeros(length(Fs_adc_range), 1);

for Fs_idx = 1:length(Fs_adc_range)
    current_Fs_adc = Fs_adc_range(Fs_idx);
    
    reconstructed_amplitudes_for_current_fs = zeros(num_events_per_amplitude, 1);
    
    for event_idx = 1:num_events_per_amplitude
        continuous_signal = representative_amplitude * H_t;
        
        % Próbkowanie ADC w losowej fazie
        max_shift_samples_analog = floor(Fs_analog / current_Fs_adc);
        random_shift_samples_analog = randi([0, max_shift_samples_analog - 1]);
        
        t_start_adc = t_analog(random_shift_samples_analog + 1);
        t_adc = t_start_adc : (1/current_Fs_adc) : T_sim;
        
        sampled_signal_clean = interp1(t_analog, continuous_signal, t_adc, 'linear', 0);
        sampled_signal_noisy = sampled_signal_clean + noise_amplitude_rms * randn(size(sampled_signal_clean));
        
        % WAŻNE: Filtr FIR powinien być zaprojektowany dla konkretnej Fs_adc.
        % Jeśli Fs_adc się zmienia, filtr 1/H(f) również powinien się zmienić!
        % Dla uproszczenia (ale z niedokładnością) na razie używamy tego samego filtru.
        % W rzeczywistości, dla każdej Fs_adc należałoby ponownie zaprojektować `impulse_response_fir`.
        
        % Poniżej, dla poprawnej analizy, MUSISZ przeliczyć filtr dla każdej Fs_adc!
        % Krok do samodzielnego wykonania:
        % 1. Dyskretne H(t) dla danej Fs_adc
        % 2. Dyskretne H(f) dla danej Fs_adc
        % 3. Projekt filtra 1/H(f) dla danej Fs_adc
        
        % Tymczasowo, dla ilustracji, używamy filtra z Fs_adc_selected.
        % To pokaże, jak ważne jest dopasowanie filtru do Fs_adc.
        
        % Aby to zrobić poprawnie, musisz ponownie wykonać kroki 2 i 6 w tej pętli,
        % dla `current_Fs_adc`.

        % PRZYKŁADOWY, POPRAWNY SPOSÓB PRZEPROWADZENIA TEGO KROKU:
        % 1. Zprobkuj H_t dla current_Fs_adc:
        H_t_sampled_for_current_Fs = interp1(t_analog, H_t, 0:(1/current_Fs_adc):T_sim, 'linear', 0);
        
        % 2. Oblicz H(f) dla tego próbkowania
        H_f_sampled_for_current_Fs = fft(H_t_sampled_for_current_Fs);
        
        % 3. Zaprojektuj nowy filtr Wienera dla current_Fs_adc
        impulse_response_fir_long_current_Fs = ifft(conj(H_f_sampled_for_current_Fs) ./ (abs(H_f_sampled_for_current_Fs).^2 + regularization_k));
        impulse_response_fir_current_Fs = fftshift(impulse_response_fir_long_current_Fs);
        impulse_response_fir_current_Fs = real(impulse_response_fir_current_Fs(floor(end/2) - floor(fir_len/2) + 1 : floor(end/2) + ceil(fir_len/2)));
        
        reconstructed_signal_raw = conv(sampled_signal_noisy, impulse_response_fir_current_Fs);
        
        delay_samples = floor(length(impulse_response_fir_current_Fs) / 2);
        reconstructed_signal_aligned = reconstructed_signal_raw(delay_samples + 1 : end - delay_samples);
        if isempty(reconstructed_signal_aligned)
            reconstructed_amplitudes_for_current_fs(event_idx) = 0;
        else
            reconstructed_amplitudes_for_current_fs(event_idx) = max(reconstructed_signal_aligned);
        end
        
    end
    mean_reconstructed_sampling(Fs_idx) = mean(reconstructed_amplitudes_for_current_fs);
    std_reconstructed_sampling(Fs_idx) = std(reconstructed_amplitudes_for_current_fs);
    reconstruction_errors_sampling(Fs_idx) = std(reconstructed_amplitudes_for_current_fs - representative_amplitude);
    
    fprintf('  Fs_ADC: %.1f MHz, Średnia zrekonstruowana: %.2f, Std dev: %.4f, Błąd: %.4f\n', ...
            current_Fs_adc / 1e6, mean_reconstructed_sampling(Fs_idx), ...
            std_reconstructed_sampling(Fs_idx), reconstruction_errors_sampling(Fs_idx));
end

%% 9. Wizualizacja wyników dla częstotliwości próbkowania
figure;
yyaxis left;
plot(Fs_adc_range/1e6, reconstruction_errors_sampling, 'o-');
ylabel('Błąd Rekonstrukcji (Std Dev)');
yyaxis right;
plot(Fs_adc_range/1e6, std_reconstructed_sampling ./ mean_reconstructed_sampling * 100, 'x-'); % Procentowa rozdzielczość energetyczna
ylabel('Procentowa Rozdzielczość [%]');
title('Błąd Rekonstrukcji vs. Częstotliwość Próbkowania ADC (dla amplitudy %.2f)', representative_amplitude);
xlabel('Częstotliwość Próbkowania ADC [MHz]');
grid on;


%% Dodatkowe wizualizacje (przykładowe)

% Przykład pojedynczego impulsu i jego rekonstrukcji
figure;
subplot(3,1,1);
plot(t_analog * 1e9, continuous_signal);
title('Oryginalny sygnał H(t)');
xlabel('Czas [ns]');
ylabel('Amplituda');
grid on;

subplot(3,1,2);
plot(t_adc * 1e9, sampled_signal_noisy);
title('Próbkowany sygnał z szumem (losowa faza)');
xlabel('Czas [ns]');
ylabel('Amplituda');
grid on;

subplot(3,1,3);
% Czas dla zrekonstruowanego sygnału
t_reconstructed = (0:length(reconstructed_signal_raw)-1) * (1/Fs_adc_selected);
plot(t_reconstructed * 1e9, reconstructed_signal_raw);
hold on;
stem((t_start_adc + (delay_samples * (1/Fs_adc_selected))) * 1e9, max(reconstructed_signal_aligned), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
title('Zrekonstruowany sygnał po dekonwolucji');
xlabel('Czas [ns]');
ylabel('Amplituda');
grid on;