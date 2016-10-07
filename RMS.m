function juice_RMS
clear all; close all;
tic;
%%%%%%%%%%%%%%%%%%%%       滤波器参数    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 250;  % Sampling Frequency    
%%% filter_design(Fs,Fstop1,Fpass1,Fpass2,Fstop2)    % 滤波器系数的计算，并存储系数 B
load('B.mat'); 
%%%%%%%%%%%%%%%%%%%%        参数       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unit_num = 50;    %多少个点算一次rms

global data_length;
data_length=[];
N2_dir_path = 'D:\spindle\RMS\Data';  % （N2_dir_path：所有被试的数据文件夹的路径）
N2_dir = dir(N2_dir_path); %（ N2_dir：N2数据文件夹的信息的结构体）
m = length(N2_dir); % m-2 为被试个数

%%%%%%%%%%%%%%%%%%%   rms, threshold   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for subj_i = 17 % 对被试的循环 
    subj_name = N2_dir(subj_i).name; % （subj_name：当前被试的名字）
    subj_dir_path = [N2_dir_path,'\',subj_name]; % （subj_dir_path：当前被试文件夹的路径） 
    subj_dir_s = dir(subj_dir_path); % (subj_dir_N23: 当前被试的文件夹的信息的结构)
    m1 =  length(subj_dir_s); % m1-2 当前被试文件加下的mat个数（即有多少段的数据）
    if isempty(data_length)
        data_length=zeros(m-2,m1-2);%建一个数组存放每个被试每个mat文件的数据长度
    else
    end
    cd(subj_dir_path);   % 
    mkdir('results');    % 在被试文件夹下建文件夹 results，用于存储结果
    subj_results_path = [subj_dir_path,'\','results'];
    
    all_mat_rms2 = [];
    
    for subj_s = 3:m1
        N_subj_mat=subj_s-2;
        subj_s_name = subj_dir_s(subj_s).name; %文件中的mat文件信息
        cd(subj_results_path);%进入results文件
%         mkdir(subj_s_name);  % 在results文件夹下建文件夹
        %cd('D:\MATLAB\matlab2012a\mywork\mywork\sunjb\RMS-A\chun_mat');%回到最初目录，否则会找不到后面要用到的函数   
        subj_results_s_path = subj_results_path;
        subj_s_path = subj_dir_path; % （subj_dir_path：当前被试文件夹的路径） 
        subj_mat = dir([subj_dir_path '\' '*.mat']);
        num_subj_mat = length(subj_mat); % （num_subj_mat：当前被试的mat个数）
        all_mat_rms = [];
        all_mat_time = zeros(1,num_subj_mat);
        all_mat_time_num = zeros(1,num_subj_mat);
        %%%%%%%%%%%%%%%%%%%%%           rms          %%%%%%%%%%%%%%%%%%%%%%%
        
        subj_mat_path = [subj_s_path '\' subj_s_name]; % （subj_mat_dir：当前被试的第N_subj_mat 个mat文件路径）
        %%%%  进入 .mat   %%%%
        load(subj_mat_path); 
        data_length(subj_i-2,subj_s-2)=length(b(3,:));
        mat_orgn_data= b([1:8,13:17],:); % （mat_orgn_data:load的最原始的数据）
        mat_orgn_length = size(mat_orgn_data,2);   % （mat_length：mat_orgn_data的长度）
        each_mat_rms = mat_rms(mat_orgn_data,unit_num,B,N_subj_mat,subj_results_s_path);
        mat_time = size(each_mat_rms,1)*0.25;   % 计算每个mat的时间,s
        all_mat_time(N_subj_mat) = mat_time;     % 保存每个mat的时间在一个变量中    
        all_mat_rms = [all_mat_rms;each_mat_rms];  % 保存每个mat的rms在一个变量中
        all_mat_time_num(N_subj_mat) = mat_orgn_length;
        all_mat_rms2 = [all_mat_rms2;all_mat_rms];
        path_2 = strcat(subj_results_s_path,'\','all_mat_time.mat');
        save(path_2,'all_mat_time');   % 保存每个stage的time，保存 all_mat_time 到 all_mat_time.mat
        save(path_2,'all_mat_time_num');
        thres= prctile(all_mat_rms2,95);  %阈值
%         clear all_mat_rms2;
        cd(subj_results_path);
        save(['thres_' subj_s_name],'thres');  % 保存 thres 到 thres.mat
        if ~exist('spindle_a0.mat')
            spd_detection_a0=zeros(data_length(subj_i-2,subj_s-2),1,m1-2,15);
        else
        end
        filter_mat_name = dir(['.\','filter_data_*.mat']);% 用来通配符 'rms_data_*.mat'
        load('all_mat_time');
        s_time = sum(all_mat_time)/60; % 该阶段的时间，min
        Spindle_general = [];
        for pole=1:11
            Spd_amplitude = [];  %保存每个电极的amplitude
            Spd_frequency = [];   % 保存每个电极的 frequency
            Spd_duration = [];   %  保存每个电极的 duration
            cd(subj_results_s_path);
            filter_data_path = [filter_mat_name(N_subj_mat).name]; % 
            load(filter_data_path);  % 导入滤波数据
            each_mat_rms_pole=each_mat_rms(:,pole);%取一个电极的数据
            rms_data_1 = (each_mat_rms_pole>thres(pole)) ;  % 把rms中大于thres的标记为1   
            res = lianxu_1(rms_data_1,unit_num,Fs);    % 计算连续1的位置和个数                     
            eval(['Sres_',num2str(N_subj_mat),'mat_',num2str(pole),'pole','=','res',';']);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%   singular spindle    %%%%%%%%%%%%%%%%%%%%%%
            cd(subj_results_s_path);   % 涉及到数据保存
            if isempty(res)
                continue;
            else
                [spd_num,unuseful] = size(res);
                for spd_i = 1:spd_num
                    st = res(spd_i,1);  % 计算rms纺锤波的起点index
                    ed =res(spd_i,1)+res(spd_i,2)-1;  % 计算rms纺锤波的结束index
                    st = (st-1)*unit_num+1;   % 换算纺锤波的起点index
                    ed = ed*unit_num;     % 换算纺锤波的结束index
                %%%%%%%%%储存纺锤波的起止点，用于画图%%%%%%%%%%%%%%%%%
                    spd_detection_a0(st:ed,1,N_subj_mat,pole)=1;
                    if ~exist('spindle_a0.mat')
                        spd_a0_position=zeros(spd_num,2,m1-2,15);                          
                        spd_a0_position(spd_i,:,N_subj_mat,pole)=[st ed];
                    else
                        spd_a0_position(spd_i,:,N_subj_mat,pole)=[st ed];
                    end
                        save('spindle_a0','spd_a0_position','spd_detection_a0');
                        spd_data = mat_data_3(st:ed);  %得到纺锤波数据 
                        [pks_p,locs_p] = findpeaks(spd_data);  % 得到极大值
                        [pks_n,locs_n] = findpeaks(-spd_data);  %得到极小值
                        locs_p = st-1+locs_p;  %得到极大值在mat_data_3 的绝对index
                        locs_n = st-1+locs_n;  %得到极小值在mat_data_3 的绝对index


           %%%%%%%%%%%%%%%     计算幅值   %%%%%%%%%%%%%%%%%%%%%%%%%       

                        length_p = length(locs_p);
                        length_n = length(locs_n);
                        length_union = length_p+length_n;
                        pks_union = zeros(1,length_union);

                        if (locs_p(1)<locs_n(1))  % 极大值在前面
                            pks_union(1:2:end) = pks_p;   
                            pks_union(2:2:end) = pks_n;    % 进行峰值的插缝拼接，便于计算pks-to-pks                       
                            pks_union_1 = [pks_union(2:end),0];   % 从第二个开始取
                            pks_union_2 = pks_union+pks_union_1;  %  计算 pks-to-pks
                            pks_union = pks_union_2(1:end-1);   % 取的实际的pks-to-pks
                            [amplitude,index_1] = max(pks_union);  % 得到最大的peak-to-peak


                        else
                            pks_union(1:2:end) = pks_n;   
                            pks_union(2:2:end) = pks_p;
                            pks_union_1 = [pks_union(2:end),0]; 
                            pks_union_2 = pks_union+pks_union_1;  %  计算 pks-to-pks
                            pks_union = pks_union_2(1:end-1); 
                            [amplitude,index_1] = max(pks_union);  % 得到最大的peak-to-peak

                        end                       
                        Spd_amplitude = [Spd_amplitude,amplitude];

            %%%%%%%%%%%%%%%%%%  计算频率  %%%%%%%%%%%%%%%%%
%                         pks_p_diff = diff(locs_p);
%                         pks_n_diff = diff(locs_n);
%                         pks_diff = [pks_p_diff pks_n_diff];
%                         frequency = Fs/mean(pks_diff);
                        N=512;
                        FFT_data=abs(fft(spd_data,N));
                        [pks,locs]=findpeaks(FFT_data(1:length(FFT_data)/2));
                        freqs=(1:N/2)*Fs/N;
                        locs_n=[];
                        freqs_n=[];
                        for i_pks=1:length(pks)
                            i_freqs=freqs(locs(i_pks));%峰值点对应的频率
                            if i_freqs>=12&&i_freqs<=16.5%找出满足频率为12-16hz的频率
                                locs_n=[locs_n locs(i_pks)];%满足频率范围的位置
                                freqs_n=[freqs_n i_freqs];
                            else
                            end
                        end     
                        if ~isempty(freqs_n)
                            frequency=max(freqs_n);%12-16hz范围内的最高频率
                            freq_floor=floor(frequency);
                            diff_freq=abs(frequency-freq_floor);
                            if(0<=diff_freq&&diff_freq<0.125)
                                frequency=freq_floor;
                            elseif(0.125<=diff_freq&&diff_freq<0.375)
                                frequency=freq_floor+0.25;
                            elseif(0.375<=diff_freq&&diff_freq<0.625)
                                frequency=freq_floor+0.5;
                            elseif(0.625<=diff_freq&&diff_freq<0.875)
                                frequency=freq_floor+0.75;
                            elseif(diff_freq>=0.875)
                                frequency=freq_floor+1;
                            end
                        else
                            for i_pks=1:length(pks)
                                i_freqs=freqs(locs(i_pks));%峰值点对应的频率
                                if i_freqs>=11&&i_freqs<=17%找出满足频率为12-16hz的频率
                                    locs_n=[locs_n locs(i_pks)];%满足频率范围的位置
                                    freqs_n=[freqs_n i_freqs];
                                else
                                end
                            end
                            frequency=max(freqs_n);%11-17hz范围内的最高频率
                            freq_floor=floor(frequency);
                            diff_freq=abs(frequency-freq_floor);
                            if(0<=diff_freq&&diff_freq<0.125)
                                frequency=freq_floor;
                            elseif(0.125<=diff_freq&&diff_freq<0.375)
                                frequency=freq_floor+0.25;
                            elseif(0.375<=diff_freq&&diff_freq<0.625)
                                frequency=freq_floor+0.5;
                            elseif(0.625<=diff_freq&&diff_freq<0.875)
                                frequency=freq_floor+0.75;
                            elseif(diff_freq>=0.875)
                                frequency=freq_floor+1;
                            end
                        end
                        Spd_frequency = [Spd_frequency,frequency];
                        
            %%%%%%%%%%%%%%%%%  计算持续时间   %%%%%%%%%%%%%%%%%%%%%%%%%
                        duration = res(spd_i,2)*(unit_num/Fs);   %  
                        Spd_duration = [Spd_duration,duration];

                end

            end

            
%             save(['Spindle_',num2str(pole),'pole'],'spd_a0_position','spd_detection_a0');           
            Spd_num = numel(Spd_duration);
            Spd_destiny= Spd_num/s_time; % 计算density=spindles/min
            mean_amplitude = mean(Spd_amplitude);
            mean_frequency = mean(Spd_frequency);
            mean_duration = mean(Spd_duration);
            general.num=Spd_num;
            general.destiny=Spd_destiny;
            general.amplitude=Spd_amplitude;
            general.frequency=Spd_frequency;
            general.duration=Spd_duration;
            general.mean_amplitude=mean_amplitude;
            general.mean_frequency=mean_frequency;
            general.mean_duration=mean_duration;
            general.data_duration=s_time;
            clear Spd_*pole;
            clear Spd_destiny;   
            clear Spd_amplitude;
            clear Spd_frequency;
            clear Spd_duration;               
            save(['Spindle_detail',subj_s_name(1:end-4),'_',num2str(pole),'pole'] ,'general');
            clear general;
        end
%         save('res_all','Sres_*');    %保存所有以'Sres_'开头的变量到 res_all.mat，其实就是res
        clear Sres_*;
        
%         delete('filter*.mat');    
        delete('rms_data*.mat');
            
    end
            
    
end
toc;

end



%%%%%%%%%%%%%%%%%%%%    filter_design    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  filter_design(Fs,Fstop1,Fpass1,Fpass2,Fstop2)

% Fs = 500;  % Sampling Frequency    
% 
% Fstop1 = 10;               % First Stopband Frequency
% Fpass1 = 10.2;             % First Passband Frequency
% Fpass2 = 16.8;            % Second Passband Frequency
% Fstop2 = 17;              % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation   %%% 把这个参数忽略成既定参数
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.001;           % Second Stopband Attenuation
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
B  = firpm(N, Fo, Ao, W, {dens});
save('B.mat','B')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%   each_mat_rms   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function each_mat_rms = mat_rms(mat_orgn_data,unit_num,B,N_subj_mat,subj_results_s_path)
mat_orgn_length = size(mat_orgn_data,2);   % （mat_length：mat_orgn_data的长度）
mat_length = mat_orgn_length-mod(mat_orgn_length,unit_num); % (mat_length:计算数据长度中能被unit_num整除的部分)
mat_data_1 = mat_orgn_data(:,1:mat_length); % 将当前数据不能被unit_num整除的部分从最后去掉 
mat_data_2 = zeros(11,mat_length); % 默认为C3电极，暂不考虑双极导联  （mat_data_2: 减去参考电压的数据 ）
mat_data_2(1,:) = mat_data_1(1,:) - mat_data_1(13,:);  %%FP1相对于A2的信号；
mat_data_2(2,:) = mat_data_1(2,:) - mat_data_1(12,:);  %%FP2相对于A1的信号；
mat_data_2(3,:) = mat_data_1(3,:) - mat_data_1(13,:);  %%F3相对于A2的信号；
mat_data_2(4,:) = mat_data_1(4,:) - mat_data_1(12,:);  %%F4相对于A1的信号；
mat_data_2(5,:) = mat_data_1(5,:) - mat_data_1(13,:);  %%C3相对于A2的信号；2
mat_data_2(6,:) = mat_data_1(6,:) - mat_data_1(12,:);  %%C4相对于A1的信号；
mat_data_2(7,:) = mat_data_1(7,:) - mat_data_1(13,:);  %%P3相对于A2的信号；
mat_data_2(8,:) = mat_data_1(8,:) - mat_data_1(12,:);  %%P4相对于A1的信号；
mat_data_2(9,:) = mat_data_1(9,:) - (mat_data_1(12,:)+mat_data_1(13,:))/2;  %%Fz相对于(A1+A2)/2的信号；
mat_data_2(10,:) = mat_data_1(10,:) - (mat_data_1(12,:)+mat_data_1(13,:))/2;  %%Cz相对于(A1+A2)/2的信号；
mat_data_2(11,:) = mat_data_1(11,:) - (mat_data_1(12,:)+mat_data_1(13,:))/2;  %%Pz相对于(A1+A2)/2的信号；
clear mat_orgn_data;
clear mat_data_1;

%%%%%%%%%%%%%%%%%%%%  滤波   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_data_2 = mat_data_2';
mat_data_3 = filtfilt(B,1,mat_data_2); % 对mat_data_2进行滤波，得到10列数据；需要较长时间，所以save  mat_data_3。

name_1 = strcat('filter_data_',num2str(N_subj_mat,'%02d'),'.mat');
path =strcat(subj_results_s_path,'\',name_1);
save(path,'mat_data_3');           % 保存滤波后的数据 mat_data_3 到 filter_data_...mat
mat_data_4 = reshape( mat_data_3,unit_num,[],11); % 对数据进行维度转换，便于一次性rms计算。
clear mat_data_2;
clear mat_data_3;

%%%%%%%%%%%%%%%%%%%%  计算rms   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
each_mat_rms= rms(mat_data_4);  %%求rms (1*10 | 1*10 | 1*10 ...)
each_mat_rms = reshape(each_mat_rms,[],11);
%each_mat_rms = each_mat_rms';
name_2 = strcat('rms_data_',num2str(N_subj_mat,'%02d'),'.mat');  % num2str(N_subj_mat,'%02d')
path =strcat(subj_results_s_path,'\',name_2);
save(path,'each_mat_rms');   % 保存each_mat_rms 到  rms_data_...mat
%%%%%%%%%%%%%%%%%%   存储当前被试的rms  %%%%%%%%%%%%%%%%%%%%%%%%%

end


%%%%%%%%%%%%%%%%%%%%%%%     lianxu_1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = lianxu_1(a,Num,Fs)
%a = [0 ones(1,5)  0 0 1 0 ones(1,10) 0] ;
b = find(a) ; %找出rms计算结果中大于阈值被赋值为1的位置
res = [] ;
n = 1 ; 
i = 1 ;
while i < length(b)
    j = i+1 ;
    while j <= length(b) &&  b(j)==b(j-1)+1 %如果出现连续为1
        n = n + 1 ;%连续为1的个数
        j = j + 1 ;
    end
    if n >= ((0.4*Fs)/Num) && n <= ((2*Fs)/Num)    % 连续为1的个数范围在此之内，则满足纺锤波的持续时间，就将这个值记录下来
        res = [res ; b(i),n] ;
   
    end
    n = 1 ;
    i = j ;
end
end
 
