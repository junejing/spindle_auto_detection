function spindle_detection_MP
warning off all
%%  参数读取
% num_subj = 30;                                                 %% 读取被试人数；
fs = 500;                         %% 采样频率；
ds = 4;                              %% 下采样比例（默认4）；
d_fs = fs/ds;                                                               %% 计算下采样后的采样率；
N = 1024;                               %% 读取每次分解的信号点数，默认
load('GaborAtom_1024.mat');                    %% 读取原子库；                                                                %% 释放不再使用的存储空间；
ep = 0.05;                                                                  %% 信号重构残差截止阈值；
H_psd = 16;                            %% 读取定义spindle的最高频率；
L_psd = 9;                             %% 读取定义spindle的最低频率；
TH1 = 7.5;                                 %% 读取定义的spindle幅度阈值（大于7.5）；
S_duration = 0.5;                           %% 读取定义spindle的最短持续时间；
L_duration = 2;                           %% 读取定义spindle的最长持续时间；
bbb = [-3.98923231507507e-05,-3.46161293505339e-05,0.00109113359144369,-0.000817126039535110,-0.00743905493062784,0.0106040711249938,0.0246986155741091,...
    -0.0577878721767207,-0.0483106137321757,0.298035851567139,0.559999006947751,0.298035851567139,-0.0483106137321757,-0.0577878721767207,0.0246986155741091,...
    0.0106040711249938,-0.00743905493062784,-0.000817126039535110,0.00109113359144369,-3.46161293505339e-05,-3.98923231507507e-05];
%%  数据读取和计算
total_len = 0;
data_path='D:\spindle\MP\data\N2';
data_dir=dir(data_path);
num_subj=length(data_dir);%N2文件夹中mat文件个数
detection=zeros(num_subj-2,300000);
for sub_matj = 4  %%对mat文件的循环
    spindle_num = zeros(1,10);
    data_len = 0;
    spindle_para = zeros(1,6,10);
    spindle_data = zeros(1,N,10);
    spindle_atom = zeros(1,N,10);
    spindle_marker = zeros(1,N,10);
    spindle_density = zeros(4,10);
    mat_path=fullfile(data_path,data_dir(sub_matj).name);
    load(mat_path);%导入mat文件
    data_matrix =b([1:10 29:32],:);  %% 读取当前数据段的有效数据；
    data_matrix = downsample(data_matrix',ds)';
    SigLen = size(data_matrix,2)-mod(size(data_matrix,2),N);  %% 计算数据长度中能被数据单元整除的部分；
    data_matrix = data_matrix(:,1:SigLen);  %% 将当前数据不能被计算单元整除的部分从最后去掉；
    data_10channels = zeros(10,SigLen);  %% 建立用于分析的数据存储空间；
    data_10channels(1,:) = data_matrix(1,:) - data_matrix(12,:);  %% FP1相对于A2的信号；
    data_10channels(2,:) = data_matrix(2,:) - data_matrix(11,:);  %% FP2相对于A1的信号；
    data_10channels(3,:) = data_matrix(3,:) - data_matrix(12,:);  %% F3相对于A2的信号；
    data_10channels(4,:) = data_matrix(4,:) - data_matrix(11,:);  %% F4相对于A1的信号；
    data_10channels(5,:) = data_matrix(5,:) - data_matrix(12,:);  %% C3相对于A2的信号；
    data_10channels(6,:) = data_matrix(6,:) - data_matrix(11,:);  %% C4相对于A1的信号；
    data_10channels(7,:) = data_matrix(7,:) - data_matrix(12,:);  %% P3相对于A2的信号；
    data_10channels(8,:) = data_matrix(8,:) - data_matrix(11,:);  %% P4相对于A1的信号；
    data_10channels(9,:) = data_matrix(9,:) - data_matrix(12,:);  %% O1相对于A2的信号；
    data_10channels(10,:) = data_matrix(10,:) - data_matrix(11,:);      %% O2相对于A1的信号；
    clear data_matrix;                                                  %% 释放不再使用的存储空间；
    num_cc = SigLen/N;
    cc = clock;
%     disp(['正在计算第',num2str(sub_matj),'个被试（共',num2str(num_subj),'个被试）的第',num2str(ll),...
%         '段数据（共',num2str(num_yuanqiang),'段数据）——数据长度',num2str(num_cc*N/d_fs),'秒——当前时间：',...
%         num2str(cc(2)),'月',num2str(cc(3)),'日—',num2str(cc(4)),'时',num2str(cc(5)),'分',num2str(round(cc(6))),'秒']); %%在MATLAB主窗口实时显示计算到了哪个被试的哪段数据；
    data_len = data_len + num_cc;           %% 初始化为：data_len = zeros(1,num_subj)；
    data_filter = filtfilt(bbb,1,data_10channels');                       %% 对数据进行35Hz的低通滤波；
    clear data_10channels;
    data_tran = reshape(data_filter,N,num_cc,[]);                  %% 将数据转换为三维：N行*num_cc列*10层
    clear data_filter;
%     spm_progress_bar('Init',num_cc*10,'Inner');                         %% 打开spm小窗口的进度条，分辨率为“1/段数*电极数”；
    guan = 0;
    for chs = 5
        for len = 1:num_cc
            guan = guan+1;
            x_xdata = data_tran(:,len,chs)';
%             TH1= prctile(x_xdata,95);
            [g_atomCoe_totle,g_atom,a,x_error] = matchPursuit(x_xdata,g,ep);
            [N1,M1] = size(g_atom);                                     %% N代表分解后选择原子的个数，M表示每个原子的长度
              for i = 1:N1
                 g_atom(i,:) = g_atom(i,:)*g_atomCoe_totle(i);
                 [atom_peak,index_atom] = findpeaks(g_atom(i,:));
                 if numel(atom_peak) < 2
                     continue;
                 end
                 yu = abs(fftshift(fft(g_atom(i,:)))).^2;               %% 周期图计算原子能量谱（正负半轴各N/2个点）；
                 yu1 = yu(N/2+1:end);                                   %% 只要正半轴的N/2个点，对应频率：0Hz 到 (N/2-1)/(N/d_fs) Hz;
                 [yu1_peak,index1] = findpeaks(yu1);                    %% 计算能量谱的峰值
                 if isempty(index1)                                     %% 检测是否有峰值，如没有峰值，计算下一个原子；
                     continue;
                 end
                 [yu1_max,index2] = max(yu1_peak);                      %% 计算能量谱的最大峰值
                 index = index1(index2);
                 index_tran0 = ((0:N-1)-N/2)*(d_fs/N);
                 index_tran = index_tran0(N/2+1:end);
                 if (index_tran(index) >= L_psd) && (index_tran(index) <= H_psd)   %% 检测能量谱的最大峰值是否在spindle频带范围内
                     B = find(g_atom(i,:) >= TH1);
                     if isempty(B)
                         continue;
                     end
                     duration = (B(end) - B(1))/d_fs;                         
                     if duration >= S_duration && duration <= L_duration
                         marker_s = zeros(1,N);
                         spindle_num(chs) = spindle_num(chs) +1;        %% 初始化为：spindle_num = zeros(10,num_subj);
                         spindle_para(spindle_num(chs),1,chs) = index_tran(index);
                         spindle_para(spindle_num(chs),2,chs) = max(atom_peak);
                         spindle_para(spindle_num(chs),3,chs) = duration;                             
                         marker_s(B(1):B(end)) = 1;
                         spd_st=(B(1)+N*(len-1)-1)*ds;
                         spd_ed=(B(end)+N*(len-1)-1)*ds;
                         detection(sub_matj-2,spd_st:spd_ed)=1;
                         spindle_data(spindle_num(chs),:,chs) = x_xdata;
                         spindle_atom(spindle_num(chs),:,chs) = g_atom(i,:);
                         spindle_marker(spindle_num(chs),:,chs) = marker_s;
                         zzz = sum(reshape(spindle_atom(spindle_num(chs),:,chs).^2,8,[]))./sum(reshape((spindle_data(spindle_num(chs),:,chs)-mean(spindle_data(spindle_num(chs),:,chs))).^2,8,[]));
                         yyy = reshape(repmat(zzz,8,1),1,[]).*spindle_marker(spindle_num(chs),:,chs);
                         duration1 = sum(yyy > 0.4)/d_fs;
                         duration2 = sum(yyy > 0.2)/d_fs;
                         spindle_para(spindle_num(chs),4,chs) = duration2; 
                         spindle_para(spindle_num(chs),5,chs) = duration1; 
                         if duration1 >= 0.3
                             spindle_para(spindle_num(chs),6,chs) = 100;
                         elseif duration1 < 0.3 && duration2 >= 0.5
                             spindle_para(spindle_num(chs),6,chs) = 50;
                         else
                             spindle_para(spindle_num(chs),6,chs) = 25;
                         end
                     end
                 end
             end              
%             spm_progress_bar('Set',guan);                               %% 对spm小窗口的进度条加1；
        end
    end
    spindle_density(1,:) = spindle_num./(data_len*N/d_fs/60);
    spindle_density(2,:) = sum(spindle_para(:,6,:) == 25)./(data_len*N/d_fs/60);
    spindle_density(3,:) = sum(spindle_para(:,6,:) == 50)./(data_len*N/d_fs/60);
    spindle_density(4,:) = sum(spindle_para(:,6,:) == 100)./(data_len*N/d_fs/60);
    save(['Sub',num2str(sub_matj-2),'_spindle_all.mat'],'spindle_num','data_len','spindle_density','spindle_para','spindle_data','spindle_atom','spindle_marker');
    total_len = total_len + data_len;
end
save('MP_detection','detection');
disp('o(∩_∩)o');
disp('o(∩_∩)o');
disp('o(∩_∩)o');
disp('o(∩_∩)o');
disp('o(∩_∩)o');
disp(['终于计算完成了！累计计算数据：',num2str(total_len*N/d_fs),'秒——当前时间：',...
    num2str(cc(2)),'月',num2str(cc(3)),'日—',num2str(cc(4)),'时',num2str(cc(5)),'分',num2str(round(cc(6))),'秒']);
% spm_progress_bar('Clear');  %% 关闭spm小窗口的进度条；







%% -----------------------------------------------------
function [g_atomCoe_totle,g_atom,a,x_error] = matchPursuit(x,g,ep)
%% 进行匹配追踪
% num = size(g,1);
num = 1000;
[row,colomn] = size(x);
if row > colomn
    x = x';
end
x_error = x;
g_atom = zeros(num,colomn);
g_atomCoe_totle = zeros(1,num);
a = zeros(1,num);
for m = 1:num
    inner_product = g*x_error';%将数据和原子做内积，取数值最大
    [inner_max_abs,index] = max(abs(inner_product));
    g_atom(m,:) = g(index,:);%每次迭代后得到的原子
    g_atomCoe_totle(m) = inner_product(index);%每次迭代的内积系数
    x_error = x_error - inner_product(index).*g(index,:);%残差
    a(m) = (norm(x_error)^2)/(norm(x)^2);
    if a(m) < ep %
        break;
    end
end
g_atom = g_atom(1:m,:);
g_atomCoe_totle = g_atomCoe_totle(1:m);
a = a(1:m);

