close all
clear
clc

%% (1) Thiết lập mô hình tín hiệu, mô hình mảng UCA 2-D và thông số mô phỏng
Fs = 48000; Ts=1/Fs; % Fsample: tần số (Hz), Tsample: chu kỳ lấy mẫu (s)
t0=Ts:Ts:0.1; % Thời gian lấy mẫu (s)
f1=1350; s0_A=1*sin(2*pi*f1*t0).'; % Signal_A (mV)
f2=1375; s0_B=1*sin(2*pi*f2*t0).'; % Signal_B (mV)
f3=1325; s0_C=1*sin(2*pi*f3*t0).'; % Signal_C (mV)
n0=length(t0); % Snapshots - NoS (Number of Samples) - Độ dài tín hiệu
SNR=30; % Tỷ số tín hiệu trên nhiễu
angles_A=-120*pi/180; % Góc tới Signal_A azimuth
angles_B=0*pi/180; % Góc tới Signal_B
angles_C=70*pi/180; % Góc tới Signal_C
c = 343; % Tốc độ lan truyền âm thanh trong không khí
lambda_max=c/min([f1 f2 f3]); % Lambda max
D=3; % Số lượng tín hiệu tới
k=2*pi/lambda_max; % Hệ số hoạt động của mảng
k1=2*pi/(c/f1); % Hệ số góc mỗi tín hiệu
k2=2*pi/(c/f2); % Hệ số góc mỗi tín hiệu
k3=2*pi/(c/f3); % Hệ số góc mỗi tín hiệu
phi=-180:0.05:180; % Góc quét phi để tạo ma trận lái mảng
t1=Ts:Ts:5; %thoi gian lay mau (s)
s1_A=1*sin(2*pi*f1*t1).'; %signal_A (mV)
s1_B=1*sin(2*pi*f2*t1).'; %signal_B (mV)
s1_C=1*sin(2*pi*f3*t1).'; %signal_C (mV)
s1=s1_A+s1_B+s1_C; %mixed signal
n1=length(t1); %bien thoi gian - snapshots
mu = 0.001; %learning rate/ toc do hoc/ he so thich nghi
Ne_0=4:2:16; % Số lượng phần tử mảng
mean_PPAR = zeros(1,length(Ne_0));
DOA_RMSE = zeros(1,length(Ne_0));
for Ne_idx = 1:length(Ne_0)
    Ne = Ne_0(Ne_idx);
    R=0.5*Ne*lambda_max/(2*pi); % Bán kính mảng
    for kk=1:1000
        a0_A=zeros(Ne,1); % Tạo ma trận 0 kích thước Nex1
        a0_B=zeros(Ne,1); % Tạo ma trận 0 kích thước Nex1
        a0_C=zeros(Ne,1); % Tạo ma trận 0 kích thước Nex1
        for l = 1:Ne
            a0_A(l)=exp(1j*k1*R*cos(angles_A-2*pi*((l-1)/Ne)));
            % Vector lái Signal_A
            a0_B(l)=exp(1j*k2*R*cos(angles_B-2*pi*((l-1)/Ne)));
            % Vector lái Signal_B
            a0_C(l)=exp(1j*k3*R*cos(angles_C-2*pi*((l-1)/Ne)));
            % Vector lái Signal_C
        end
        x0_A = zeros(Ne,n0); % Tạo ma trận 0 kích thước Nexn0
        x0_B = zeros(Ne,n0); % Tạo ma trận 0 kích thước Nexn0
        x0_C = zeros(Ne,n0); % Tạo ma trận 0 kích thước Nexn0
        for l = 1:Ne
            x0_A(l,:)=a0_A(l)*s0_A; % x0_A=a0_A*s0_A tín hiệu A tại mảng
            x0_B(l,:)=a0_B(l)*s0_B; % x0_B=a0_B*s0_B tín hiệu B tại mảng
            x0_C(l,:)=a0_C(l)*s0_C; % x0_C=a0_C*s0_C tín hiệu C tại mảng
        end
        x0=x0_A+x0_B+x0_C; %x0=∑x0i (i=1,2,...,Ne); Tổng tín hiệu tại mảng
        x0=awgn(x0,SNR,'measured'); % Thêm nhiễu Gaussian trắng vào tập tín hiệu
        %% (2) Ap dung thuat toan MUSIC
        Rx=x0*x0'/n0; % Tạo ma trận hiệp phương sai
        [eigvec,eigval]=eig(Rx); % Tính giá trị riêng và vector riêng
        En=eigvec(:,1:Ne-D); % Xây dựng ma trận vector riêng của nhiễu
        a0=zeros(Ne,1); % Tạo ma trận 0 kích thước Nex1
        Pmusic=zeros(1,length(phi));
        % Tạo ma trận 0 kích thước length(phi)
        for pp=1:length(phi) %quet phi tu -180:180
            for l=1:Ne
                a0(l,1)=exp(1j*k*R*cos(phi(pp)*pi/180-2*pi*((l-1)/Ne)));
                %tao ma tran vector lai cua mang
            end
            Pmusic(1,pp)=abs(1/(a0'*En*(En')*a0)); %cong thuc thuat toan MUSIC
        end
        Pmusic1 = imregionalmax(Pmusic);
        [biendoPmusic,index]=maxk(Pmusic(Pmusic1),D); %xac dinh goc toi tin hieu
        vitrigoctoi = find(Pmusic1);
        vitrigoctoi = vitrigoctoi(index);
        DOA_est = sort(phi(vitrigoctoi)); %qua thuat toan MUSIC
        angles_A1=DOA_est(1); %xac dinh goc toi tin hieu A
        angles_B1=DOA_est(2); %xac dinh goc toi tin hieu B
        angles_C1=DOA_est(3); %xac dinh goc toi tin hieu C
        PPAR(kk) = sum(biendoPmusic)/mean(Pmusic);
        DOA_error_A = angles_A*180/pi - angles_A1;
        DOA_error_B = angles_B*180/pi - angles_B1;
        DOA_error_C = angles_C*180/pi - angles_C1;
        RMSE(kk) = sqrt(DOA_error_A^2 + DOA_error_B^2 + DOA_error_C^2);
    end
    mean_PPAR(1,Ne_idx) = 10*log10(mean(PPAR));
    DOA_RMSE(1,Ne_idx) = mean(RMSE);
end