A = 1;%amplituda
tau = 1e-6;%stała czasowa [s]
n = 1;%rząd shapera

t_max = 10 * tau;
dt = 1e-8;
t = 0:dt:t_max;

%Odpowiedź  
H = (t/tau).^n .* exp(-t/tau);
H = H / max(H);

figure; plot(t, H); xlabel('Czas [s]'); ylabel('H(t)'); title('Odpowiedź impulsowa shapera');

%Sygnał na wyjściu
Vobs = A * H;

%Dodanie szumu 
noise_amp = 0.01;
Vobs_noisy = Vobs + noise_amp * randn(size(Vobs));

figure; plot(t, Vobs_noisy);
xlabel('Czas [s]'); ylabel('V_{obs}(t)'); title('Sygnał obserwowany z szumem');

%Próbkowanie ADC z losowym offsetem 
fs = 10e6; Ts = 1/fs;
t0 = rand * Ts;              
ts = t0:Ts:t_max;

V_samples  = interp1(t, Vobs_noisy, ts, 'pchip', 0);% próbki z szumem
V_samples0 = interp1(t, A*H       , ts, 'pchip', 0);% próbki wzorcowe

figure; plot(t, Vobs_noisy); hold on;
stem(ts, V_samples, 'r');
xlabel('Czas [s]'); ylabel('Sygnał');
title('Próbkowanie z losowym offsetem fazowym');
legend('V_{obs}(t)', 'V[n]');

%Filtr odwrotny projektowany na SIATCE ADC
ts0  = 0:Ts:t_max;
h_adc = interp1(t, H, ts0, 'pchip', 0);

%Długości i FFT
N  = numel(V_samples);
Lh = numel(h_adc);
M  = 2^nextpow2(N + Lh - 1);

H_d = fft(h_adc, M);
lambda = 1e-4;%regularyzacja
G = conj(H_d) ./ (abs(H_d).^2 + lambda);%filtr Wiener/Tikhonov

V_samples  = V_samples(:).';
V_samples0 = V_samples0(:).';
h_adc      = h_adc(:).';

%Długości i FFT
N  = numel(V_samples);
Lh = numel(h_adc);
M  = 2^nextpow2(N + Lh - 1);%zero-padding

H_d = fft(h_adc, M);

%Filtr odwrotny Wiener/Tikhonov
lambda = 1e-3;
G = conj(H_d) ./ (abs(H_d).^2 + lambda);

%Dekonwolucja przez mnożenie w czestotliwości (z paddingiem)
Xn = fft([V_samples,  zeros(1, M-N)]);
X0 = fft([V_samples0, zeros(1, M-N)]);

Y  = Xn .* G;
Y0 = X0 .* G;

V_rec  = real(ifft(Y));   V_rec  = V_rec(1:N);%wynik dla danych zaszumionych
V_rec0 = real(ifft(Y0));  V_rec0 = V_rec0(1:N);%tor "idealny" do kalibracji

[amp_est, k] = max(V_rec);
if k>1 && k<length(V_rec)
    y1=V_rec(k-1); y2=V_rec(k); y3=V_rec(k+1);
    denom = (y1 - 2*y2 + y3);
    if abs(denom)>eps
        d = 0.5*(y1 - y3)/denom;%przesunięcie sub-próbkowe
        amp_est = y2 - 0.25*(y1 - y3)*d;%amplituda w wierzchołku
    else
        d = 0;
    end
else
    d = 0;
end
t_est = ts(k) + d*Ts;

fprintf('A=%.3f, A_hat=%.3f, błąd = %.2f %% , t_hat = %.3f us\n', ...
        A, amp_est, 100*(amp_est/A-1), t_est*1e6);

%Kalibracja skali i korekcja baseline'u
scale = A / max(V_rec0);
V_rec = scale * V_rec;
bwin  = max(1, N-20):N;
V_rec = V_rec - median(V_rec(bwin));

%Estymacja amplitudy z poprawką paraboliczną
[amp_est, k] = max(V_rec);
if k>1 && k<length(V_rec)
    y1=V_rec(k-1); y2=V_rec(k); y3=V_rec(k+1);
    denom = (y1 - 2*y2 + y3);
    if abs(denom)>eps
        d = 0.5*(y1 - y3)/denom;
        amp_est = y2 - 0.25*(y1 - y3)*d;
    else
        d = 0;
    end
else
    d = 0;
end
t_est = ts(k) + d*Ts;

figure; stem(V_rec,'filled');
xlabel('n'); ylabel('d_{rec}[n]');
title(sprintf('Dekonwolucja w dziedzinie częstotliwości (Wiener), \\lambda=%.1e', lambda));
fprintf('Oszacowana amplituda: %.4f (A=%.4f), t_est=%.3f us\n', amp_est, A, t_est*1e6);
figure;
f = (0:M-1)/M * fs;
semilogx(f, 20*log10(abs(H_d)+1e-12), 'LineWidth',1.2); hold on;
semilogx(f, 20*log10(abs(G)+1e-12),   'LineWidth',1.2);
semilogx(f, 20*log10(abs(H_d.*G)+1e-12), 'LineWidth',1.2);
grid on; xlabel('f [Hz]'); ylabel('amplituda [dB]');
legend('|H|','|G|','|H\cdot G| \approx 1'); 
title(sprintf('Charakterystyki: shaper i filtr odwrotny (\\lambda=%.1e)', lambda));

lams = logspace(-4,-1,20);
abs_err = zeros(size(lams));
for i = 1:numel(lams)
    Gi = conj(H_d) ./ (abs(H_d).^2 + lams(i));
    V0i = real(ifft( fft([V_samples0 zeros(1,M-N)]).*Gi )); V0i = V0i(1:N);
    Vi  = real(ifft( fft([V_samples  zeros(1,M-N)]).*Gi )); Vi  = Vi(1:N);
    Vi  = (A/max(V0i)) * Vi;
    Vi  = Vi - median(Vi(max(1,N-20):N));
    [Ai, k] = max(Vi);
    if k>1 && k<length(Vi)
        y1=Vi(k-1); y2=Vi(k); y3=Vi(k+1);
        denom = (y1 - 2*y2 + y3); 
        if abs(denom)>eps
            d = 0.5*(y1 - y3)/denom; 
            Ai = y2 - 0.25*(y1 - y3)*d;
        end
    end
    abs_err(i) = abs(Ai - A);
end
[~,idx] = min(abs_err);
fprintf('Optymalna lambda (single shot): %.3e\n', lams(idx));
figure; loglog(lams, abs_err, '-o'); grid on;
xlabel('\lambda'); ylabel('|A_{hat}-A|'); title('Dobór regularyzacji');