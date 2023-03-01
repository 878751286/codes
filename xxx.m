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
N = 1024; % 子载波数
N_data = N; %信息序列长度
N_cp = 500;%循环前缀CP的长度
fs = 8e3;%采样率
M=2;
fl = 20e3;
fh = 30e3;
fc = 1.5e3; %载波中心频率

B = N/T_OFDM ;%信号带宽为10kHz

t = 0:1/fs:T_OFDM-1/fs;	% 一个符号周期的时间矢量

%% 常规方法生成OFDM信号
% 复基带信号
k = -N/2:1:N/2-1;
for kk = 1:1:N
    fsig(kk) = fc+k(kk)/T_OFDM;
end

% 载波信号
carrier  = zeros(N,length(t));
for k = 1:N
	carrier(k,:)  = exp(1j*2*pi*fsig(k)*t);
end

%生成信息序列
data = randi([0,1],1,N_data);%K0个子载波上要发射的信号
data = 2*data -1; % 变成1或者-1，虽然不知道为啥

% 滚降滤波器参数
alpha = 0.7;
span = 6; % 截断符号范围 
hrc = rcosdesign(alpha,span,1); % 升余弦滚降滤波器
hrc_order = length(hrc)-1; % 滚降滤波器阶数
x_bb = conv(hrc,data);  % 成型滤波
x_bb = x_bb(hrc_order/2+1:end-hrc_order/2);

x_pb = zeros(1,length(t));
for k = 1:1:N
	x_pb = x_pb+x_bb(k).*carrier(k,:);% 子载波叠加
end

x_pb = 2*real(x_pb);

figure(10);
plot(t,x_pb);
xlabel('Time(s)');
ylabel('Amp.');
title('时域信号');

nfft = 102400;
X_pb = fft(x_pb,nfft)*2/length(x_pb);
f = (0:1:nfft/2-1)*fs/nfft;
figure(20);
plot(f,abs(X_pb(1:nfft/2)));
xlabel('频率(Hz)');
ylabel('Amp.');
title('频域');



