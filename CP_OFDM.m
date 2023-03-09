%#################################################
% 程序功能：生成CP-OFDM信号
% 创建人：wangshan
% 创建时间：2023/03/07
%#################################################
clc;
clear all;
close all;
%% =================基本参数设置=================
Rb = 1e3;%比特率
fs = 8e3;% 采样频率
Ts = 1/fs;%采样间隔
fc = 2e3;%载波调制频率
B = Rb;%带宽就是最大信息率，bps
N_data = 1024;%让子载波数目=有效数据序列长度
T_data = N_data/Rb;%持续时间为1.024s
len_data = T_data*fs;
N_cp = 128;%循环前缀的长度

%生成信息序列
data = randi([0 1],1,N_data);
data = 2*data-1;

% 滚降滤波器参数
beta = 0.7; % 滚降系数
span = 6; % 截断的符号数
sps = 1; % 一个符号(这里是一个bit)采一个点
hrc = rcosdesign(beta,span,sps); % 滚降滤波器系数
hrc_order = length(hrc)-1;%滚降滤波器阶数

%% ==================基带调制===========================
% 加窗,升余弦滚降窗
data_win= conv(hrc,data); 
data_win = data_win(hrc_order/2+1:end-hrc_order/2);

% 一个OFDM符号
t1 = 0:Ts:T_data-1/fs; % 一个OFDM符号的时间轴
ofdm = zeros(1,len_data);
for k= 1:N_data %第k个子载波
    ofdm_k = data_win(k).*exp(1j*2*pi*k/T_data*t1); %把信息序列调制到子载波上
    ofdm = ofdm_k+ofdm; %各分路时域的叠加 
end % end of for k

% 加循环前缀  组帧[cp ofdm ofdm]
cp = ofdm(end-N_cp/Rb*fs+1:end);
cp_ofdm = [cp ofdm ofdm];
%% =====================通带调制乘载波============================
t1 = (0:length(cp_ofdm)-1)/fs; 
xt = 2*real(exp(1j.*2.*pi.*fc.*t1).*cp_ofdm);
figure(01);
subplot(2,1,1);
plot(t1,xt); %绘制通带信号波形，ms为单位
title('时域波形');xlabel('\bf Time(ms)');ylabel('\bf Amplitude');

% 对信号频谱分析
Nfft = 2^nextpow2(length(xt));
Xt = fft(xt,Nfft)/length(xt);
Xt = Xt(1:Nfft/2);%结果的前半部分对应正频率部分
freq_axis = (0:Nfft/2-1)/Nfft*fs; % 频率轴(Hz)
subplot(2,1,2);
plot(freq_axis,abs(Xt));
title('频谱分析');xlabel('\bf Freq(Hz)');ylabel('\bf Amplitude(V)');

%%  ==================  IFFT方法  ================================
ifft_ofdm = ifft(data_win,N_data);% IFFT变换
ifft_cp = ifft_ofdm(end-N_cp+1:end);% 插入循环前缀
ifft_cp_ofdm =[ifft_cp ifft_ofdm ifft_ofdm]; % 组帧

% 乘载波，进行通带调制
t2 = (0:length(ifft_cp_ofdm)-1)/fs;% 时间轴按照点数来
ifft_pb_xt = 2*real(ifft_cp_ofdm.*exp(1j*2*pi*fc.*t2));%通带信号取实部
% 绘制时域波形
figure(02);
subplot(2,1,1);plot(t2,ifft_pb_xt);
title('IFFT-时域波形');xlabel('\bf Time(ms)');ylabel('\bf Amplitude');
% 进行频谱分析
ifft_Nfft = 2^nextpow2(length(ifft_pb_xt));
ifft_Xt = fft(ifft_pb_xt,ifft_Nfft)/length(ifft_pb_xt);
freq_axis =(0:ifft_Nfft/2-1)/ifft_Nfft*fs ;%频率轴
subplot(2,1,2);plot(freq_axis,abs(ifft_Xt(1:ifft_Nfft/2)));
title('IFFT-频谱分析');xlabel('\bf Freq(Hz)');ylabel('\bf Amplitude(V)');


