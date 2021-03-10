clear all; close all; clc;
pkg load symbolic


set(0, 'DefaultTextInterpreter', 'latex')

% pro syntézu signálu
pocet_harmonickych = 10;

%% Definice signálu -- obdélníkový signál, sudá funkce
h = 5;
T_0 = 4; % [s] Základní perioda signálu
tau = 0.5 * T_0; % [s] 
f_0 = 1 / T_0; % [Hz] Základní frekvence
%omega_0 = 2*pi*f_0; % [rad/s] Základní úhlová frekvence

% Vytvoøení symbolických promìnných
syms t f real;
syms k integer; % Promìnná "k" je celé èíslo

% Vytvoøení signálu 
x = h * (heaviside(t - (-tau/2)) - ...
    heaviside(t - (+tau/2)));

%%% Zobrazení signálu v èase
fig_obdelnik_symb = figure('Name', 'Signal x(t) - symbolicky'); ezplot(t, x, [-T_0/2 T_0/2]); 
xlabel('t [s]'); ylabel('x(t)'); title('x(t)');


%% Koeficienty Fourierovy øady v komplexní­m tvaru 
%c_0 = 1/T_0 * int(x, t, -T_0/2, T_0/2);
c_0 = double(1/T_0 * int(h, t, -tau/2, tau/2));
c_k = 1/T_0 * ...
    int(x * exp(-1i * 2*pi*f_0 * k * t),...
    t, -T_0/2, T_0/2);
%%% Zobrazení­ spektra signálu
c_k_num = ... % Numerické hodnoty koeficientù c_k
    double([subs(c_k, k, [-pocet_harmonickych:-1]), ...
    c_0,...
    subs(c_k, k, [1:pocet_harmonickych])]);
vektor_kf_0 = f_0*[-pocet_harmonickych:pocet_harmonickych];

fig_spektrum_c_k_symb = figure('Name', 'Spektrum c[f_0 k] - symbolicky'); 
subplot(2, 1, 1); 
stem(vektor_kf_0, abs(c_k_num)); axis([-f_0*pocet_harmonickych f_0*pocet_harmonickych 0 h/tau]); hold on;
% Prùbìh spektra, pokud by "k" bylo realné èí­slo 
kk = -pocet_harmonickych : 1/10 : pocet_harmonickych;
c_kk_num = ...
    double([subs(c_k, k, kk(kk<0)), ...
    c_0,...
    subs(c_k, k, kk(kk>0))]);
plot(kk*f_0, abs(c_kk_num), 'r--');
xlabel('f [Hz]'); ylabel('|c[f_0 k]|');
XTickLabel_String = [];

subplot(2, 1, 2); 
stem(vektor_kf_0, angle(c_k_num)); axis([-f_0*pocet_harmonickych f_0*pocet_harmonickych -3.5 3.5]);
xlabel('f [Hz]'); ylabel('c[f_0 k]');

fig_spektrum_c_k = figure('Name', 'Spektrum c[f_0 k]'); 
subplot(2, 1, 1); 
stem(vektor_kf_0, abs(c_k_num)); axis([-f_0*pocet_harmonickych f_0*pocet_harmonickych 0 h/tau]); hold on;
plot(kk*f_0, abs(c_kk_num), 'r--');
xlabel('f [Hz]'); ylabel('|c[f_0 k]|');
subplot(2, 1, 2); 
stem(vektor_kf_0, angle(c_k_num)); axis([-f_0*pocet_harmonickych f_0*pocet_harmonickych -3.5 3.5]);
xlabel('f [Hz]'); ylabel('c[f_0 k]');


print(fig_spektrum_c_k, '-dpdf', 'obdelnikovy_signal_FR_spektrum_c_k.pdf');
print(fig_spektrum_c_k, '-dsvg', 'obdelnikovy_signal_FR_spektrum_c_k.svg');

rc_kk_num = real(c_kk_num); rc_kk_num(abs(rc_kk_num) < 1e-10) = 0;
ic_kk_num = imag(c_kk_num); ic_kk_num(abs(ic_kk_num) < 1e-10) = 0;
c_kk_num = rc_kk_num + 1i*ic_kk_num; clear rc_kk_num ic_kk_num;

% Uložení­ promìnných do souboru
csvwrite('obdelnikovy_signal_FR_abs_spojita_obalka_spektra.csv', 3*[f_0*kk.' abs(c_kk_num).']);
csvwrite('obdelnikovy_signal_FR_abs_koeficienty_c_k.csv', 3*[vektor_kf_0.' abs(c_k_num).']);
csvwrite('obdelnikovy_signal_FR_arg_koeficienty_c_k.csv', [3*vektor_kf_0.' angle(c_k_num).']);

%%% Syntéza signálu z koeficientù Fourierovy øady v komplexní­m tvaru
fig_synteza_c = figure('Name', 'Synteza z c[f_0 k]');
x_synt_kompl = c_0; % Stejnosmìrná složka
for kk = 1:pocet_harmonickych 
    x_synt_kompl = x_synt_kompl + ...
        subs(c_k, k, -kk) * ... % conj(subs(c_k, k, kk)) * ...
        exp(1i*2*pi*f_0*(-kk)*t) + ...
        subs(c_k, k, kk) * ...
        exp(1i*2*pi*f_0*kk*t)
    ezplot(t, x_synt_kompl, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}(t)');
    title(['Pocet harmonickych = ', num2str(kk)]);
end

print(fig_synteza_c, '-dpdf', 'obdelnikovy_signal_synteza_c_k.pdf');
print(fig_synteza_c, '-dsvg', 'obdelnikovy_signal_synteza_c_k.svg');

clear c_k_num kk vektor_kf_0 x_synt_kompl

%% Koeficienty Fourierovy øady v reálném tvaru, c_k => a_k, b_k
a_0 = 2 * c_0;
a_k =  2 * real(c_k);
b_k = -2 * imag(c_k); 
a_k_num = double(subs(a_k, k, 1:pocet_harmonickych));
b_k_num = double(subs(b_k, k, 1:pocet_harmonickych));
%%% Zobrazení spektra signálu
fig_spektrum_a_k_b_k_symb = figure('Name', 'Spektrum a[f_0 k] a b[f_0 k] - symbolicky'); 
kk = 1/10 : 1/10 : pocet_harmonickych;
aa_k_num = double([a_0,...
    subs(a_k, k, kk)]);
subplot(1, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [a_0 a_k_num]); hold on;
plot([0 kk*f_0], aa_k_num, '--r');
title('Realne spektrum');
xlabel('f [Hz]'); ylabel('a[f_0 k]');
kk = 1 : 1/10 : pocet_harmonickych;
bb_k_num = double(subs(b_k, k, kk));
subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN b_k_num]); hold on;
plot(kk*f_0, zeros(size(bb_k_num)), '--r');
xlabel('f [Hz]'); ylabel('b[f_0 k]');

fig_spektrum_a_k_b_k = figure('Name', 'Spektrum a[f_0 k] a b[f_0 k]');
kk = 1/10 : 1/10 : pocet_harmonickych;
subplot(2, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [a_0 a_k_num]); hold on;
plot([0 kk*f_0], aa_k_num, '--r');
title('Realne spektrum');
xlabel('f [Hz]'); ylabel('a[f_0 k]');

kk = 1 : 1/100 : pocet_harmonickych;
bb_k_num = subs(b_k, k, kk);
title('magnitudove spektrum');
subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN b_k_num]); hold on;
plot(kk*f_0, zeros(size(bb_k_num)), '--r');
xlabel('f [Hz]'); ylabel('b[f_0 k]');
title('magnitudove spektrum');
print(fig_spektrum_a_k_b_k, '-dpdf', 'obdelnikovy_signal_FR_spektrum_a_k_b_k.pdf');
print(fig_spektrum_a_k_b_k, '-dsvg', 'obdelnikovy_signal_FR_spektrum_a_k_b_k.svg');

%%% Syntéza signálu z koeficientù Fourierovy øady v reálném tvaru
fig_synteza_a_k_b_k = figure('Name', 'Synteza z a[f_0 k] a b[f_0 k]');
x_synt_real = a_0/2;
for kk = 1:pocet_harmonickych
    x_synt_real = x_synt_real + ...
        a_k_num(kk) * cos(2*pi*f_0*kk*t) + ...
        b_k_num(kk) * sin(2*pi*f_0*kk*t);
    ezplot(t, x_synt_real, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}(t)');
    title(['Poèet harmonických = ', num2str(pocet_harmonickych)]);
    %pause
end


print(fig_synteza_a_k_b_k, '-dpdf', 'obdelnikovy_signal_synteza_a_k_b_k.pdf');
print(fig_synteza_a_k_b_k, '-dsvg', 'obdelnikovy_signal_synteza_a_k_b_k.svg');

clear a_0 a_k a_k_num b_k b_k_num kk x_synt_real

%% Koeficienty Fourierovy øady v I. elektrotechnickém tvaru, c_k => A_k, phi_k
A_0 = c_0; % Stejnosmìrná složka 
A_k = 2 * abs(c_k); % Amplitudy kosinusovek
A_k_num = subs(A_k, k, 1:pocet_harmonickych);
phi_k_num = -double(angle( subs(c_k, k, 1:pocet_harmonickych) )); % Fázový posuv kosinusovek
%%% Zobrazení­ spektra signálu
fig_spektrum_A_k_phi_k_symb = figure('Name', 'Spektrum A[f_0 k] a phi[f_0 k] - symbolicky'); 
kk = 1/10 : 1/10 : pocet_harmonickych;
AA_k_num = [2*A_0,...
    double(subs(A_k, k, kk))];
subplot(2, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [A_0 double(A_k_num)]); hold on;
plot([0 kk*f_0], AA_k_num, '--r');
title('Analytické spektrum');
xlabel('f [Hz]'); ylabel('A[f_0 k]');


subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN phi_k_num]); hold on;
xlabel('f [Hz]'); ylabel('phi[f_0 k]');


fig_spektrum_A_k_phi_k = figure('Name', 'Spektrum A[f_0 k] a phi[f_0 k]');
subplot(2, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [A_0 double(A_k_num)]); hold on;
plot([0 kk*f_0], AA_k_num, 'r--');
title('Analytické spektrum');
xlabel('f [Hz]'); ylabel('A[f_0 k]');
subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN phi_k_num]);
xlabel('f [Hz]'); ylabel('phi[f_0 k]');

print(fig_spektrum_A_k_phi_k, '-dpdf', 'obdelnikovy_signal_FR_spektrum_A_k_phi_k.pdf');
print(fig_spektrum_A_k_phi_k, '-dsvg', 'obdelnikovy_signal_FR_spektrum_A_k_phi_k.svg');

% Syntéza signálu z koeficientù Fourierovy øady v I. elektrotechnickém tvaru
fig_synt_A_phi = figure('Name', 'Synteza z A[f_0 k] a phi[f_0 k]');
x_synt_I_elt = double(A_0);
for kk = 1:pocet_harmonickych
    if(isnan(phi_k_num(kk)))
  x_synt_I_elt = x_synt_I_elt + ...
        double(A_k_num(kk)) * ...
        cos(2*pi*f_0 * kk *t);
  else
  x_synt_I_elt = x_synt_I_elt + ...
        double(A_k_num(kk)) * ...
        cos(2*pi*f_0 * kk *t - phi_k_num(kk));
   endif
    ezplot(t, x_synt_I_elt, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}(t)');
    title(['Poèet harmonických = ', num2str(pocet_harmonickych)]);
    %pause
end


print(fig_synt_A_phi, '-dpdf', 'obdelnikovy_signal_synteza_A_k_phi_k.pdf');
print(fig_synt_A_phi, '-dsvg', 'obdelnikovy_signal_synteza_A_k_phi_k.svg');

clear A_0 A_k A_k_num kk phi_k_num x_synt_I_elt

%return

%%%% Koeficienty Fourierovy øady v II. elektrotechnickém tvaru, c_k => A_k, psi_k
A_0 = c_0; % Stejnosìrná složka
A_k = 2 * abs(c_k); % Amplitudy kosinusovek
A_k_num = subs(A_k, k, 1:pocet_harmonickych);
psi_k_num = -double(angle( subs(1i * c_k, k, 1:pocet_harmonickych) )); % Fázový posuv sinusovek
kk = 1/10 : 1/10 : pocet_harmonickych;
AA_k_num = [A_0,...
    double(subs(A_k, k, kk))];
%%%% Zobrazení spektra signálu
figure;
subplot(2, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [A_0 double(A_k_num)]); hold on;
plot([0 kk*f_0], AA_k_num, 'r--');
title('Reálné spektrum');
xlabel('f [Hz]'); ylabel('A[k]');
subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN psi_k_num]);
xlabel('f [Hz]'); ylabel('psi[k]');
%%%% Syntéza signálu z koeficientù Fourierovy øady v I. elektrotechnickém tvaru
figure;
x_synt_II_elt = A_0;
for kk = 1:pocet_harmonickych
  if(isnan(psi_k_num(kk)))
  x_synt_II_elt = x_synt_II_elt + ...
        A_k_num(kk) * ...
        sin(2*pi*f_0 * kk *t);
  else
    x_synt_II_elt = x_synt_II_elt + ...
        A_k_num(kk) * ...
        sin(2*pi*f_0 * kk *t - psi_k_num(kk));
   endif
    ezplot(t, x_synt_II_elt, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}');
    title(['Poèet harmonických = ', num2str(pocet_harmonickych)]);
    %pause
end
clear A_0 A_k A_k_num kk psi_k_num x_synt_II_elt
clear c_0 c_k

%%%% Koeficienty Fourierovy øady v reálném tvaru
a_0 = 2 * double(1/T_0 * int(x, t, -tau/2, tau/2));
a_k = 2 * 1/T_0 * ...
    int(x * cos(2*pi*f_0 * k *t), ...
    t, -tau/2, tau/2);
b_k = 2 * 1/T_0 * ...
    int(x * sin(2*pi*f_0 * k *t), ...
    t, -tau/2, tau/2);
%%%% Zobrazení spektra signálu
figure;
subplot(2, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [a_0 double(subs(a_k, k, 1:pocet_harmonickych))]);
title('Reálné spektrum');
xlabel('f [Hz]'); ylabel('a[k]');
subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN double(subs(b_k, k, 1:pocet_harmonickych))]);
xlabel('f [Hz]'); ylabel('b[k]');
%%%% Syntéza signálu z koeficientù Fourierovy øady v reálném tvaru
figure;%14
x_synt_real = a_0 / 2;
for kk = 1:pocet_harmonickych
    x_synt_real = x_synt_real + ...
      subs(a_k, k, kk) * cos(2*pi*f_0 * kk * t) + ...
      subs(b_k, k, kk) * sin(2*pi*f_0 * kk * t)
    ezplot(t, x_synt_real, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}');
    title(['Poèet harmonických = ', num2str(pocet_harmonickych)]);
    %pause
end
clear kk x_synt_real

%%%% Koeficienty Fourierovy øady v komplexní­m tvaru, a_k, b_k => c_k 
c_0 = a_0 / 2;
c_k = (a_k - 1i*b_k) / 2;
%%%% Zobrazení­ spektra signálu
figure;
subplot(2, 1, 1);
stem(f_0 * [-pocet_harmonickych:pocet_harmonickych], ...
    abs([double(subs(c_k, k, [-pocet_harmonickych:-1])), ...
    c_0, ...
    double(subs(c_k, k, [1:pocet_harmonickych]))]));
title('Amplitudové spektrum');
xlabel('f [Hz]'); ylabel('|c[kf_0]|');
subplot(2, 1, 2);
stem(f_0 * [-pocet_harmonickych:pocet_harmonickych], ...
    angle([double(subs(c_k, k, [-pocet_harmonickych:-1])), ...
    c_0, ...
    double(subs(c_k, k, [1:pocet_harmonickych]))]));
title('Fázové spektrum');
xlabel('f [Hz]'); ylabel('arg(c[kf_0])');
%%%% Syntéza signálu z koeficientù Fourierovy øady v komplexním tvaru
figure;%16
x_synt_kompl = c_0;
for kk = 1:pocet_harmonickych
    x_synt_kompl = x_synt_kompl + ...
        subs(c_k, k, -kk) *...
        exp(1i * 2*pi*f_0 * (-kk) * t) + ...
        subs(c_k, k, kk) *...
        exp(1i * 2*pi*f_0 * kk * t)
    ezplot(t, x_synt_kompl, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}');
    title(['Poèet harmonických = ', num2str(pocet_harmonickych)]);
    %pause
end
clear c_0 c_k kk x_synt_kompl x_synt_real

%%%% Koeficienty Fourierovy øady v I. elektrotechnickém tvaru, a_k, b_k => A_k, phi_k
A_0 = a_0 / 2; % Stejnosmìrná složka 
A_k = sqrt(a_k^2 + b_k^2); % Amplitudy kosinusovek
A_k_num = subs(A_k, k, 1:pocet_harmonickych);
phi_k_num = -double(angle( subs((a_k - 1i*b_k) / 2, k, 1:pocet_harmonickych)) ); % Fázový posuv kosinusovek
%%%% ZobrazenÃí spektra signálu
figure;
subplot(2, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [A_0 double(A_k_num)]);
title('Reálné spektrum');
xlabel('f [Hz]'); ylabel('A[k]');
subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN phi_k_num]);
xlabel('f [Hz]'); ylabel('phi[k]');
%%%% Syntéza signálu z koeficientù Fourierovy øady v I. elektrotechnickém tvaru
figure;
x_synt_I_elt = A_0;
for kk = 1:pocet_harmonickych
  if(isnan(phi_k_num(kk)))
  x_synt_I_elt = x_synt_I_elt + ...
        A_k_num(kk) * ...
        cos(2*pi*f_0 * kk *t)
  else
    x_synt_I_elt = x_synt_I_elt + ...
        A_k_num(kk) * ...
        cos(2*pi*f_0 * kk *t - phi_k_num(kk))
  endif
    ezplot(t, x_synt_I_elt, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}');
    title(['Poèet harmonických = ', num2str(pocet_harmonickych)]);
    %pause
end
clear A_0 A_k A_k_num kk phi_k_num x_synt_I_elt

%%%% Koeficienty Fourierovy øady v II. elektrotechnickém tvaru (sin), c_k => A_k, psi_k
A_0 = a_0 / 2; % Stejnosmìrná složka 
A_k = sqrt(a_k^2 + b_k^2); % Amplitudy kosinusovek
A_k_num = subs(A_k, k, 1:pocet_harmonickych);
psi_k_num = double(angle( subs(b_k - 1i*a_k, k, 1:pocet_harmonickych)) ); % Fázový posuv sinusovek
%%%% Zobrazení spektra signálu
figure;
subplot(2, 1, 1);
stem(f_0 * [0:pocet_harmonickych], ...
    [A_0 double(A_k_num)]);
title('Reálné spektrum');
xlabel('f [Hz]'); ylabel('A[k]');
subplot(2, 1, 2);
stem(f_0 * [0:pocet_harmonickych], ...
    [NaN psi_k_num]);
xlabel('f [Hz]'); ylabel('psi[k]');
%%%% Syntéza signálu z koeficientù Fourierovy øady v I. elektrotechnickém tvaru
figure;
x_synt_II_elt = A_0;
for kk = 1:pocet_harmonickych
  if(isnan(psi_k_num(kk)))
  x_synt_II_elt = x_synt_II_elt + ...
        double(A_k_num(kk)) * ...
        sin(2*pi*f_0 * kk *t)
  else
    x_synt_II_elt = x_synt_II_elt + ...
        double(A_k_num(kk)) * ...
        sin(2*pi*f_0 * kk *t - psi_k_num(kk))
   endif
    ezplot(t, x_synt_II_elt, [-3*T_0/2 3*T_0/2]);
    xlabel('t [s]'); ylabel('x_{synt}');
    title(['Poèet harmonických = ', num2str(pocet_harmonickych)]);
    %pause
end
clear A_0 A_k A_k_num kk psi_k_num x_synt_II_elt
clear a_0 a_k b_k