% close all
% clear
% clc

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
Ne=6; % Số lượng phần tử mảng
D=3; % Số lượng tín hiệu tới
R=0.045; % Bán kính mảng
k=2*pi/lambda_max; % Hệ số hoạt động của mảng
k1=2*pi/(c/f1); % Hệ số góc mỗi tín hiệu
k2=2*pi/(c/f2); % Hệ số góc mỗi tín hiệu
k3=2*pi/(c/f3); % Hệ số góc mỗi tín hiệu
phi=-180:0.05:180; % Góc quét phi để tạo ma trận lái mảng
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
% Tạo ma trận 0 kích thước 1length(phi)
for pp=1:length(phi) %quet phi tu -180:180
    for l=1:Ne
        a0(l,1)=exp(1j*k*R*cos(phi(pp)*pi/180-2*pi*((l-1)/Ne)));
        %tao ma tran vector lai cua mang
    end
    Pmusic(1,pp)=abs(1/(a0'*En*(En')*a0)); %cong thuc thuat toan MUSIC
end

%% (3) Thiet lap thong so dau vao truoc khi ap dung Beamforming
[biendoPmusic,vitrigoctoi]=findpeaks(Pmusic); %xac dinh goc toi tin hieu
DOA_est = phi(vitrigoctoi); %qua thuat toan MUSIC
angles_A1=DOA_est(1)*pi/180; %xac dinh goc toi tin hieu A
angles_B1=DOA_est(2)*pi/180; %xac dinh goc toi tin hieu B
angles_C1=DOA_est(3)*pi/180; %xac dinh goc toi tin hieu C
t1=Ts:Ts:5; %thoi gian lay mau (s)
s1_A=1*sin(2*pi*f1*t1).'; %signal_A (mV)
s1_B=1*sin(2*pi*f2*t1).'; %signal_B (mV)
s1_C=1*sin(2*pi*f3*t1).'; %signal_C (mV)
s1=s1_A+s1_B+s1_C; %mixed signal
[Pxx, ff1] = pwelch(s1, [], [], [], Fs);
n1=length(t1); %bien thoi gian - snapshots
a1_A=zeros(Ne,1); %tao ma tran 0 kich thuoc Nex1
a1_B=zeros(Ne,1); %tao ma tran 0 kich thuoc Nex1
a1_C=zeros(Ne,1); %tao ma tran 0 kich thuoc Nex1
for l = 1:Ne
    a1_A(l)=exp(1j*k1*R*cos(angles_A1-2*pi*((l-1)/Ne))); %vector lai signal_A
    a1_B(l)=exp(1j*k2*R*cos(angles_B1-2*pi*((l-1)/Ne))); %vector lai signal_B
    a1_C(l)=exp(1j*k3*R*cos(angles_C1-2*pi*((l-1)/Ne))); %vector lai signal_C
end
x1_A=zeros(Ne,n1); %tao ma tran 0 kich thuoc Nexn0
x1_B=zeros(Ne,n1); %tao ma tran 0 kich thuoc Nexn0
x1_C=zeros(Ne,n1); %tao ma tran 0 kich thuoc Nexn0
for l = 1:Ne
    x1_A(l,:)=a1_A(l)*s1_A; %x1_A=a1_A*s1_A tin hieu A tai mang
    x1_B(l,:)=a1_B(l)*s1_B; %x1_B=a1_B*s1_B tin hieu B tai mang
    x1_C(l,:)=a1_C(l)*s1_C; %x1_C=a1_C*s1_C tin hieu C tai mang
end
x1=x1_A+x1_B+x1_C; %x0=∑x0i (i=1,2,...,Ne); tong tin hieu thu tai mang
x1=awgn(x1,SNR,'measured'); %them nhieu Gaussian trang vao tin hieu
%% (4) Ap dung thuat toan LMS
mu_0 = [0.0001 0.0005 0.001 0.005 0.01]; %learning rate/ toc do hoc/ he so thich nghi
for i = 1:length(mu_0)
    mu=mu_0(i);
    y = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
    e = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
    w = zeros(1, Ne); %tao ma tran 0 kich thuoc 1xn1
    abswn = zeros(Ne, n1); %tao ma tran 0 kich thuoc Nexn1
    wn = zeros(Ne, n1); %tao ma tran 0 kich thuoc Nexn1
    for m=1:n1
        y(m)=w*x1(:,m); %dau ra bo loc y
        e(m)=s1_C(m)'-y(m); %sai so/ loi
        w = w + mu*x1(:,m)'*e(m); %cap nhat trong so w
        abswn(:,m)=abs(w); %gia tri abs trong so w
        wn(:,m)=w; % gia tri trong so w
    end
    a1=zeros(Ne,1);
    LMS=zeros(1,length(phi));
    for the=1:length(phi) %quet phi tu -180:180
        for l=1:Ne
            a1(l,1)=exp(1j*k*R*cos(phi(the)*pi/180-2*pi*((l-1)/Ne)));
            %tao ma tran vector lai cua mang
        end
        LMS(the)=abs(w*a1); %cong thuc pho dap ung thuat toan LMS
    end
    %Bieu dien thuat toan MUSIC
    % figure (1);
    % plot(phi,20*log10(LMS/max(LMS))); hold on
    figure(2);
    plot(t1,e); hold on
    % figure(3);
    % plot(t1,abswn(1)); hold on
end
% figure(1); grid on; legend('mu=0.0001','mu=0.0005','mu=0.001','mu=0.005','mu=0.01');
figure(2); grid on; legend('mu=0.0001','mu=0.0005','mu=0.001','mu=0.005','mu=0.01');