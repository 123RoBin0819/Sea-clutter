abs_data_matrix = abs(amplitude_complex_T1);

%时频图和B显图
%{
power_matrix = abs_data_matrix.^2; % 计算功率
% 转换为分贝
db_matrix = 10 * log10(power_matrix);
% 绘制图像
figure;
imagesc(db_matrix);
colormap('summer'); % 设置颜色映射，这里使用'jet' colormap
colorbar; % 添加颜色条
title('T1脉冲回波'); % 添加标题
xlabel('采样点数'); % x轴标签
ylabel('脉冲数'); % y轴标签

% 选择第445列数据
data_col_100 = abs_data_matrix(:, 100);
data_col_445 = abs_data_matrix(:, 445);
data_col_616 = abs_data_matrix(:, 616);

% 设置STFT参数
window_length = 256;
overlap = 200;
fs = 1; % 虚拟单位频率，因为只进行频率分析

% 计算短时傅里叶变换
[s, f, t] = spectrogram(data_col_100, hamming(window_length), overlap, window_length, fs);

% 将正频率结果镜像得到负频率结果
s_neg = conj(flipud(s(2:end,:))); % 镜像正频率

% 合并正负频率结果
s_combined = [s_neg; s];
f_combined = [-flipud(f(2:end)); f];

% 绘制时频谱图（包括正负频率）
figure;
subplot(2, 1, 1);
surf(t/2000, f_combined, 10*log10(abs(s_combined)), 'EdgeColor', 'none');
view(0, 90); % 设置视角
xlabel('时间/s');
ylabel('频率');
title('Short-Time Fourier Transform (STFT) 目标位置：100');
colormap('jet'); % 设置颜色映射，这里使用'jet' colormap
colorbar; % 添加颜色条

% 设置STFT参数
window_length = 256;
overlap = 200;
fs = 1; % 虚拟单位频率，因为只进行频率分析

% 计算短时傅里叶变换
[s, f, t] = spectrogram(data_col_445, hamming(window_length), overlap, window_length, fs);

% 将正频率结果镜像得到负频率结果
s_neg = conj(flipud(s(2:end,:))); % 镜像正频率

% 合并正负频率结果
s_combined = [s_neg; s];
f_combined = [-flipud(f(2:end)); f];

% 绘制时频谱图（包括正负频率）
subplot(2, 1, 2);
surf(t/2000, f_combined, 10*log10(abs(s_combined)), 'EdgeColor', 'none');
view(0, 90); % 设置视角
xlabel('时间/s');
ylabel('频率');
title('Short-Time Fourier Transform (STFT) 目标位置：445');
colormap('jet'); % 设置颜色映射，这里使用'jet' colormap
colorbar; % 添加颜色条

% 设置STFT参数
window_length = 256;
overlap = 200;
fs = 1; % 虚拟单位频率，因为只进行频率分析

% 计算短时傅里叶变换
[s, f, t] = spectrogram(data_col_616, hamming(window_length), overlap, window_length, fs);

% 将正频率结果镜像得到负频率结果
s_neg = conj(flipud(s(2:end,:))); % 镜像正频率

% 合并正负频率结果
s_combined = [s_neg; s];
f_combined = [-flipud(f(2:end)); f];




% 绘制图像
data_col_100 = abs_data_matrix(1:2048, 100);
data_col_445 = abs_data_matrix(1:2048, 445);
t=1:2048;
%绘制二维图
figure;
subplot(2, 1, 1);
plot(t/2000,data_col_100);
title('采集点100，纯海杂波');
xlabel('时间/s');
ylabel('幅度');

subplot(2, 1, 2);
plot(t/2000,data_col_445);
title('采集点445，含目标波');
xlabel('时间/s');
ylabel('幅度');

%}








%总体幅度图
%{
pulse_points = 1:size(abs_data_matrix, 1);
columns = 1:size(abs_data_matrix, 2);
% 定义每个数据段的长度
segment_length = 1024;

% 计算数据段的数量
num_segments = floor(length(pulse_points) / segment_length);

% 初始化平均值数组
average_values = zeros(num_segments, size(abs_data_matrix, 2));

% 计算每个数据段的平均值
for i = 1:num_segments
    start_idx = (i - 1) * segment_length + 1;
    end_idx = i * segment_length;
    segment_data = abs_data_matrix(start_idx:end_idx, :);
    average_values(i, :) = mean(segment_data);
end

% 可视化平均值
figure;
imagesc(average_values);
colormap('jet'); % 设置颜色映射
colorbar; % 添加颜色条
title('以1024间隔，平均值B显图');
xlabel('Column Index');
ylabel('Segment Index');

% 创建行和列的索引
pulse_points = 1:size(average_values, 1);
columns = 20:size(average_values, 2);

% 创建网格
[PulsePoints, Columns] = meshgrid(pulse_points, columns);

% 转置维度
average_values = average_values(:,20:950)';

% 绘制三维图像
figure;
mesh(Columns, PulsePoints, average_values);
xlabel('距离采样点');
ylabel('时间间隔，以1024分割');
zlabel('幅度');
title('杂波总览');
colorbar;  % 添加颜色条
%}








%信杂比
%{
data_col_100 = abs_data_matrix(:, 200);%海杂波
data_col_445 = abs_data_matrix(:, 612);
data_col_616 = abs_data_matrix(:, 616);
signal_power = mean(data_col_445.^2, 'all');

% 计算噪声功率
noise_power = mean(data_col_100.^2, 'all');

% 计算平均信噪比（以分贝为单位）
snr_db = 10 * log10((signal_power-noise_power) / noise_power);
disp(['Average SNR (dB): ', num2str(snr_db)]);
%}








%非广延熵检测

% 指定数组大小
row = 950;
col = 512;
Fs = 2000;       % 采样频率
 
T = 1/Fs;       % 采样时间
 
L = 256;        % 信号长度
N=256;
q=2;
t = (0:L-1)*T; % 时间
s=zeros(row,col);
for i = 1:row
    for j = 1:col
        data_col = abs_data_matrix(256*(j-1)+1:256*j, i);
        Y = fft(data_col,N)/N*2;   %除以N乘以2才是真实幅值，N越大，幅值精度越高
        Y(1)=Y(1)/2;
        A = abs(Y);     %幅值
        A=A.^2/sum(A.^2);
        s(i,j)=(1-sum(A.^q))/(q-1);
    end
end
% 创建行和列的索引
row = 1:row;
col = 1:col;

% 创建网格
[row, col] = meshgrid(col, row);
surf(row, col, s, 'EdgeColor', 'none');
view(0, 90); % 设置视角
xlabel('时间单元');
ylabel('位置单元');
title('非广延熵比较 q=0.5');
colormap('jet'); % 设置颜色映射，这里使用'jet' colormap
colorbar; % 添加颜色条

n_cs=0;
n_ts=0;

for j=1:512
    threshold = func_CACFAR(s(:,j),j);
    for i=1:950
        if((i<608||i>619)&&s(i,1)<threshold(i,1))
            n_cs=n_cs+1;
        elseif((i>=608&&i<=619)&&s(i,1)<threshold(i,1))
            n_ts=n_ts+1;
        end
    end
end
p_d=n_ts/(512*12)
p_fa=n_cs/(512*938)











%单点的非广延熵
%{
%快速傅里叶变换
data_col_449 = abs_data_matrix(1:256, 459);
Fs = 2000;       % 采样频率
 
T = 1/Fs;       % 采样时间
 
L = 256;        % 信号长度
 
t = (0:L-1)*T; % 时间
 

figure;
 
plot(t,data_col_449)
 
title('加噪声的信号')
 
xlabel('时间(s)')
 
N = 2^nextpow2(L); %采样点数，采样点数越大，分辨的频率越精确，N>=L，超出的部分信号补为0
Y = fft(data_col_449,N)/N*2;   %除以N乘以2才是真实幅值，N越大，幅值精度越高
Y(1)=Y(1)/2;
f = Fs/N*(0:1:N-1); %频率
 
A = abs(Y);     %幅值
 
P = angle(Y);   %相值
 
figure;
 
plot(f(1:N),A(1:N));   %函数fft返回值的数据结构具有对称性,因此我们只取前一半
 
title('幅值频谱')
 
xlabel('频率(Hz)')
 
ylabel('幅值')

%非广延熵

A=A.^2/sum(A.^2);
q=2;
s=(1-sum(A.^q))/(q-1);
%}







%多普勒功率谱
%{
Fs = 2000;       % 采样频率
T = 1/Fs;       % 采样时间
N = 1024;

f = Fs/N*(-N/2:1:N/2); % 频率轴
data_col_449 = amplitude_complex_T1(1:1024, 449);
b = zeros(N, 1);
S = zeros(1, N);
for j = 1:N
    Z = data_col_449.* exp(-2*pi*(1:N)*T*f(j));
    b(j,1)=sum(Z(j,:));
    S(j) = 1/N * abs(b(j,1))^2;
end
figure;
plot(f(1:N), 10*log10(S(1:N)));   
title('含目标单元');
xlabel('频率(Hz)');
ylabel('多普勒功率谱');
%}








%时间自相关
%{
m=128;
data=amplitude_complex_T1(1:m,200);
c=sum(abs(data));
acf=zeros(m,1);
for i=1:m
    b=zeros(m,1);
    for j=1:m-i
        b(j,1)=data(j,1)*conj(data(j+i-1));
    end
    acf(i,1)=abs(sum(b))/c;
end
acf=acf./acf(1,1)
i=1:m
stem(i/2,acf,'fill');
title('纯海杂波时间自相关图');
xlabel('持续时间/(ms)');
ylabel('时间相关系数');
%}





%空间自相关
%{
data=mean(abs(amplitude_complex_T1));
m=mean(data);
M=950
spacf=zeros(M,1);
for i=1:M
    b=zeros(M,1);
    for j=1:M-i
        b(j,1)=(data(1,j)-m)*(data(1,j+i-1)-m);
    end
    spacf(i,1)=sum(b);
end
spacf=spacf./spacf(1,1);
i=1:M
stem(i,spacf,'fill');
title('空间自相关性图');
xlabel('采样距离');
ylabel('空间相关系数');
%}





%平均相对幅度
%{
data=abs(amplitude_complex_T1(:,449));%目标波
cankao=abs(amplitude_complex_T1(:,300));%海杂波
zabo=abs(amplitude_complex_T1(:,200));
raa=zeros(1024,1);
for i=1:1024
    raa(i,1)=mean(data(128*(i-1)+1:128*i,1))/mean(zabo(128*(i-1)+1:128*i,1));
end
i=1:1024;
plot(i,raa);
hold on; % 保持图形，以便画下一个折线图
raa1=zeros(1024,1);
for i=1:1024
    raa1(i,1)=mean(cankao(128*(i-1)+1:128*i,1))/mean(zabo(128*(i-1)+1:128*i,1));
end
i=1:1024
plot(i,raa1);
legend('目标波', '海杂波', 'Location', 'northeast'); % 添加图例，设定位置为左上角
title('目标波相对海杂波的相对平均幅度');
xlabel('数据段（以128分隔）');
ylabel('相对平均幅度');
hold off;
%}


