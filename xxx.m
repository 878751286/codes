%#################################################
% 程序功能：生成CP-OFDM信号
% 创建人：wangshan
% 创建时间：2023/02/28
% 存在的问题：写的不对，还得改，明天再改
%#################################################
 clc;
 clear;
 close all;

%% =================基本参数设置=================
T_OFDM = 1.024; % 单个OFDM码元时间 s
num_OFDM_symbol = 1; % OFDM符号数
N = 1024; % 子载波数
N_data = num_OFDM_symbol*N; % 原始比特数
N_cp = 500; % 循环前缀CP的长度
fs = 8e3; % 采样率
M = 2; % 调试阶数
fc = 1.5e3; % 载波中心频率
B = N/T_OFDM ; % 信号带宽
sps = fs/B;

%生成信息序列
% infobits = randi([0,1],1,N_data);%K0个子载波上要发射的信号
% save ./mat_data/infobits.mat infobits;
load ./mat_data/infobits.mat infobits;
data = infobits;
data = 2*data-1; % 变成1或者-1，虽然不知道为啥

alpha = 0.7;
span = 6; % 截断符号范围 
hrc = rcosdesign(alpha,span,1); % 升余弦滚降滤波器
hrc_order = length(hrc)-1; % 滚降滤波器阶数


%% 常规方法生成OFDM信号
% 复基带信号
fsig = zeros(1,N);
k = -N/2:1:N/2-1;
for kk = 1:1:N
    fsig(kk) = fc+k(kk)/T_OFDM;
end

% 载波信号
t1 = 0:1/fs:T_OFDM-1/fs;	% 一个符号周期的时间矢量
carrier  = zeros(N,length(t1));
for k = 1:1:N
	carrier(k,:) = exp(1i*2*pi*fsig(k)*t1);
end

% 加窗
x_bb = conv(hrc,data); 
x_bb = x_bb(hrc_order/2+1:end-hrc_order/2);

x_pb = zeros(1,length(t1));
for k = 1:1:N
	x_pb = x_pb+x_bb(k).*carrier(k,:);% 子载波叠加
end

x_pb = 2*real(x_pb);

figure(10);
subplot(221);
plot(t1,x_pb);
xlabel('Time(s)');
ylabel('Amp.');
title('时域信号（定义）');

figure(10);
subplot(222);
nfft = 2^nextpow2(length(x_pb));
X_pb = fft(x_pb,nfft)*2/length(x_pb);
f1 = (0:1:nfft/2-1)*fs/nfft;
plot(f1,abs(X_pb(1:nfft/2)));
xlabel('频率(Hz)');
ylabel('Amp.');
title('频域（定义）');

%% DFT实现
NFFT = N; % FFT点数，等于子载波数
ifftSignal = ifft([data(1:NFFT/2) zeros(1,(sps-1)*NFFT,1) data(end-NFFT/2+1:end)]); % 频域插值对应时域升采样，将每个OFDM码元长度提高sps倍

t2 = (0:1:length(ifftSignal)-1)/fs;
x_pb_DFT = 2*real(ifftSignal.*exp(-1i*2*pi*fc*t2)); % 发射的应为实信号
x_pb_DFT = conv(hrc,x_pb_DFT); % 加窗
x_pb_DFT = x_pb_DFT(hrc_order/2+1:end-hrc_order/2);

figure(10);
subplot(223);
plot(t2,x_pb_DFT);
xlabel('Time(s)');
ylabel('Amp.');
title('时域信号（DFT实现）');

nfft = 2^nextpow2(length(x_pb_DFT));
X_pb_DFT = fft(x_pb_DFT,nfft)*2/length(x_pb_DFT);
f2 = (0:1:nfft/2-1)*fs/nfft;

figure(10);
subplot(224);
plot(f2,abs(X_pb_DFT(1:nfft/2)));
xlabel('频率(Hz)');
ylabel('Amp.');
title('频域（DFT实现）');




