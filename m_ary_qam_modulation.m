% Clearing previous data
clc;
clear all;
close all;

% Define M-ary for QAM modulation
M = 64;

% Input validation
fprintf('\n\n\n');

% Check if log2(M) is an integer
Ld = log2(M);
ds = ceil(Ld);
dif = ds - Ld;
if (dif ~= 0)
    error('the value of M is only acceptable if log2(M) is an integer');
end

% Generate binary information
nbit = 1024; % number of information bits
msg = round(rand(nbit, 1)); % generate binary data
disp('Binary information at transmitter:');
disp(msg);
fprintf('\n\n');

% Digital Signal Representation of Binary Information
x = msg;
bp = 0.000001; % bit period
bit = [];
for n = 1:1:length(x)
    if x(n) == 1
        se = ones(1, 100);
    else
        se = zeros(1, 100);
    end
    bit = [bit se];
end

t1 = bp / 100:bp / 100:100 * length(x) * (bp / 100);
figure(1)
subplot(3, 1, 1);
plot(t1, bit, 'lineWidth', 2.5);
grid on;
axis([0 bp * length(x) -.5 1.5]);
ylabel('Amplitude (volt)');
xlabel('Time (sec)');
title('Transmitting information as digital signal');

% Reshape binary information for M-ary QAM modulation
msg_reshape = reshape(msg, log2(M), nbit / log2(M))';
disp('Information reshaped for conversion to symbolic form:');
disp(msg_reshape);
fprintf('\n\n');
size(msg_reshape);
for (j = 1:1:nbit / log2(M))
    for (i = 1:1:log2(M))
        a(j, i) = num2str(msg_reshape(j, i));
    end
end
as = bin2dec(a);
ass = as';
figure(1)
subplot(3, 1, 2);
stem(ass, 'Linewidth', 2.0);
title('Serial symbol for M-ary QAM modulation at transmitter');
xlabel('n (discrete time)');
ylabel('Magnitude');
disp('Symbolic form information for M-ary QAM:');
disp(ass);
fprintf('\n\n');

% Mapping for M-ary QAM modulation
x1 = [0:M-1];
p = qammod(ass, M); % constellation design for M-ary QAM according to symbol
sym = 0:1:M-1; % considered symbol of M-ary QAM, just for scatter plot
pp = qammod(sym, M); % constellation diagram for M-ary QAM
scatterplot(pp), grid on;
title('Constellation diagram for M-ary QAM');

% M-ary QAM Modulation
RR = real(p);
II = imag(p);
sp = bp * 2; % symbol period for M-ary QAM
sr = 1 / sp; % symbol rate
f = sr * 2;
t = sp / 100:sp / 100:sp;
ss = length(t);
m = [];
for (k = 1:1:length(RR))
    yr = RR(k) * cos(2 * pi * f * t); % in-phase or real component
    yim = II(k) * sin(2 * pi * f * t); % quadrature or imaginary component
    y = yr + yim;
    m = [m y];
end
tt = sp / 100:sp / 100:sp * length(RR);
figure(1);
subplot(3, 1, 3);
plot(tt, m);
title('Waveform for M-ary QAM modulation according to symbolic information');
xlabel('Time (sec)');
ylabel('Amplitude (volt)');

% M-ary QAM Demodulation
m1 = [];
m2 = [];
for n = ss:ss:length(m)
    t = sp / 100:sp / 100:sp;
    y1 = cos(2 * pi * f * t); % in-phase component
    y2 = sin(2 * pi * f * t); % quadrature component
    mm1 = y1 .* m((n - (ss - 1)):n);
    mm2 = y2 .* m((n - (ss - 1)):n);
    z1 = trapz(t, mm1); % integration
    z2 = trapz(t, mm2); % integration
    zz1 = round(2 * z1 / sp);
    zz2 = round(2 * z2 / sp);
    m1 = [m1 zz1];
    m2 = [m2 zz2];
end

% Demapping for M-ary QAM modulation
clear i;
clear j;
for (k = 1:1:length(m1))
    gt(k) = m1(k) + j * m2(k);
end
gt
ax = qamdemod(gt, M);
figure(3);
subplot(2, 1, 1);
stem(ax, 'linewidth', 2);
title('Re-obtained symbol after M-ary QAM demodulation');
xlabel('n (discrete time)');
ylabel('Magnitude');
disp('Re-obtained symbol after M-ary QAM demodulation:');
disp(ax);
fprintf('\n\n');

% Representation of Receiving Binary Information as Digital Signal
x = ax;
bp = 0.000001; % bit period
bit = [];
for n = 1:1:length(x)
    if x(n) == 1
        se = ones(1, 100);
    else
        se = zeros(1, 100);
    end
    bit = [bit se];
end
t1 = bp / 100:bp / 100:100 * length(x) * (bp / 100);
figure(3)
subplot(2, 1, 2);
plot(t1, bit, 'lineWidth', 2.5);
grid on;
axis([0 bp * length(x) -.5 1.5]);
ylabel('Amplitude (volt)');
xlabel('Time (sec)');
title('Receiving information as digital signal after M-ary QAM demodulation');
