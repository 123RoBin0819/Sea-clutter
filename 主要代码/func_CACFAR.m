% 函数定义：计算恒虚警率（Constant False Alarm Rate, CFAR）阈值，并绘制信号及CFAR阈值图
function y = func_CACFAR(signal,a)
    % 定义参考窗口长度
    refLength = 30;
    
    % 定义保护间隔长度
    guardLength = 13;
    
    % 定义偏移量，通过调整偏移量改变虚警率
    offset = -0.07;%HH0.041 VV0.034
    
    % 创建滑动窗函数，中心部分为参考窗口，两侧为保护间隔
    cfarWin = ones((refLength + guardLength) * 2 + 1, 1);
    cfarWin(refLength + 1 : refLength + 1 + 2 * guardLength) = 0; % 设置保护间隔为0
    cfarWin = cfarWin / sum(cfarWin); % 归一化窗函数
    
    % 计算噪声水平，通过对信号进行加权平均（卷积操作）
    noiseLevel = conv(signal, cfarWin, 'same');
    
    % 计算CFAR阈值，基于噪声水平加上偏移量
    cfarThreshold = noiseLevel + offset;
    cfarThreshold = fillmissing(cfarThreshold, 'linear');
    % 绘制图形
    if(a==1)
        figure
        plot(signal); % 绘制原始信号
        hold on; % 保持当前图像，在同一张图上添加新的数据
        plot(cfarThreshold, 'r', 'LineWidth', 2); % 绘制CFAR阈值线，红色，线宽为2
        legend('Signal', 'CFAR Threshold'); % 图例显示
        xlabel('分辨率单元索引'); % x轴标签
        ylabel('熵值'); % y轴标签
        title('恒虚警检测阈值'); % 图形标题，显示SNR值
    end

    y=cfarThreshold;
end