% 加载MAT文件中的雷达信号数据
%load('20221112150043_stare_HH.mat'); % 替换为您的MAT文件路径

% 获取雷达信号数据
radar_signal = amplitude_complex_T1; % 假设雷达信号存储在名为'amplitude_complex_T1'的变量中
[num_samples, num_cells] = size(radar_signal); % 获取样本数量和单元数量



%小波变换
%{
% 创建存储去噪后信号的变量
radar_signal_deNoise = zeros(num_samples, num_cells);
% 小波去噪参数设置
level = 4; % 小波变换的层数
wavelet_name = 'db3'; % 小波基函数名称
threshold_type = 's'; % 阈值类型
% 对每个单元的信号进行小波去噪处理
for i = 1:num_cells
    % 获取当前单元的信号
    current_signal = radar_signal(:, i);
    
    % 计算当前信号的幅度
    amplitude = abs(current_signal);
    
    % 进行小波变换
    [c, l] = wavedec(amplitude, level, wavelet_name);
    
    % 计算阈值
    sigma = median(abs(c))/0.6745;
    threshold = sigma*sqrt(2*log(length(amplitude)));
    
    % 应用软阈值处理
    c_deNoise = wthresh(c, threshold_type, threshold);
    
    % 重构去噪后的信号
    amplitude_deNoise = waverec(c_deNoise, l, wavelet_name);
    
    % 将去噪后的幅度数据转换为复数形式，并与原始信号相乘得到去噪后的复数信号
    radar_signal_deNoise(:, i) = current_signal ./ amplitude .* amplitude_deNoise;
end
abs_data_matrix=abs(radar_signal_deNoise);
%}





%中值滤波
%{
radar_signal=abs(radar_signal);
filtered_data_complex_T1 = zeros(size(amplitude_complex_T1));
for i = 1:size(radar_signal, 2)
    signal = radar_signal(:, i);
    filtered_signal = medfilt1(signal, 3); % 这里的 3 是中值滤波的窗口大小，可以根据需要调整
    filtered_data_complex_T1(:, i) = filtered_signal;
end
abs_data_matrix=filtered_data_complex_T1;
%}






%{
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
plot(data_col_445);
title('采集点445，含目标波');
xlabel('时间/s');
ylabel('幅度');
%}





%绘制三维图
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
data_col_445 = abs_data_matrix(:, 449);
signal_power = mean(data_col_445.^2, 'all');

% 计算噪声功率
noise_power = mean(data_col_100.^2, 'all');

% 计算平均信噪比（以分贝为单位）
snr_db = 10 * log10((signal_power-noise_power) / noise_power);
disp(['Average SNR (dB): ', num2str(snr_db)]);
%}



%非广延熵计算

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
figure;
surf(row, col, s, 'EdgeColor', 'none');
view(0, 90); % 设置视角
xlabel('时间单元');
ylabel('位置单元');
title('非广延熵比较 q=2');
colormap('jet'); % 设置颜色映射，这里使用'jet' colormap
colorbar; % 添加颜色条



n_cs=0;
n_ts=0;

for j=1:512
    threshold = func_CACFAR(s(:,j),j);
    for i=1:950
        if((i<608||i>619)&&s(i,1)<threshold(i,1))%VV:608——619 HH：442——454
            n_cs=n_cs+1;
        elseif((i>=608&&i<=619)&&s(i,1)<threshold(i,1))
            n_ts=n_ts+1;
        end
    end
end
p_d=n_ts/(512*12)
p_fa=n_cs/(512*938)
