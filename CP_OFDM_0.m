%#################################################
% 程序功能：生成CP-OFDM信号
% 创建人：wangshan
% 创建时间：2023/02/28
% 存在的问题：写的不对，还得改，明天再改
%#################################################
  clear all;
  close all;
  clc;

%% =================基本参数设置=================
N_data = 1000; %信息序列长度
N_cp = 500;%循环前缀CP的长度
fs = 100e3;%采样率
Rb = 5e3;%比特率
M=2;
Rs = Rb/log2(M);
T0 = N_data/Rs;%信息序列持续时间,为0.2s
Tcp = N_cp/Rs;%循环前缀CP持续时间,为0.1s
fl = 20e3;
fh = 30e3;
fc = 25e3; %载波中心频率
B =fh-fl ;%信号带宽为10kHz
K0 = B*T0;%子载波数为2e3

%生成信息序列
data = randi([0,1],1,N_data);%K0个子载波上要发射的信号
data = 2*data -1; % 变成1或者-1，虽然不知道为啥
%生成CP前缀
cp = randi([0,1],1,N_cp);
cp = 2*cp-1;
%生成后面补的0
N_zeros = 500;
T_zeros = N_zeros/Rb;%补的0序列持续时间
add0 = zeros(1,N_zeros);

% 低通滤波器参数
lpf_order = 256;%低通滤波器阶数
lpf_wn = 2*30e3/fs; %截止频率为30kHz, 
lpf_bb = fir1(lpf_order, lpf_wn);%设计好的低通滤波器

%% 基带CP-OFDM信号
data_all = [data data];%需要传输的所有信息 为什么不是单独给每一个ofdm符号用一个for循环最后再组合在一起
for k=1:K0
    trans_data(k) = data_all(k)*exp(sqrt(-1)*2*pi*k/(2*T0));
end % end of for k 
signal = [cp trans_data ];%组帧

% 经过低通滤波
signal = conv(signal,lpf_bb);
signal = signal(lpf_order/2+1:end-lpf_order/2);

%% 通带CP-OFDM信号
t = (0:length(signal)-1)/fs*1e3;
xt = 2*real(signal.*exp(sqrt(-1).*2.*pi.*fc.*t));
figure(02);
plot(t,xt); %绘制通带信号波形
xlabel('\bf Time(ms)');
ylabel('\bf Amplitude');

%对信号频谱分析
Nfft = 2^nextpow2(length(xt));
Xt = fft(xt,Nfft)/length(xt);
Xt = Xt(1:Nfft/2);
figure(03);
freq_axis = (0:Nfft/2-1)/Nfft*fs/1e3; % 频率轴(kHz)
plot(freq_axis,abs(Xt));

