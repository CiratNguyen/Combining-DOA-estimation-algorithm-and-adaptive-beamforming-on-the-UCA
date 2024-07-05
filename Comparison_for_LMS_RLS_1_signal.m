% close all
% clear
% clc

%% (1) Thiết lập mô hình tín hiệu, mô hình mảng UCA 2-D và thông số mô phỏng
Fs = 48000; Ts=1/Fs; % Fsample: tần số (Hz), Tsample: chu kỳ lấy mẫu (s)
t1=Ts:Ts:5; % Thời gian lấy mẫu (s)
f1=1350; s0_A=1*sin(2*pi*f1*t0).'; % Signal_A (mV)
SNR=30; % Tỷ số tín hiệu trên nhiễu
angles_A=0*pi/180; % Góc tới Signal_A azimuth
c = 343; % Tốc độ lan truyền âm thanh trong không khí
lambda_max=c/f1; % Lambda max
Ne=6; % Số lượng phần tử mảng
D=1; % Số lượng tín hiệu tới
R=0.5*Ne*lambda_max/(2*pi); % Bán kính mảng
k=2*pi/lambda_max; % Hệ số hoạt động của mảng
k1=2*pi/(c/f1); % Hệ số góc mỗi tín hiệu
phi=-180:0.05:180; % Góc quét phi để tạo ma trận lái mảng
% t1=Ts:Ts:5; %thoi gian lay mau (s)
s1_A=1*sin(2*pi*f1*t1).'; %signal_A (mV)
s1=s1_A;
n1=length(t1); %bien thoi gian - snapshots
mu = 0.001; %learning rate/ toc do hoc/ he so thich nghi
a1_A=zeros(Ne,1); %tao ma tran 0 kich thuoc Nex11
for l = 1:Ne
    a1_A(l)=exp(1j*k1*R*cos(0-2*pi*((l-1)/Ne))); %vector lai signal_A
end
x1_A=zeros(Ne,n1); %tao ma tran 0 kich thuoc Nexn0
for l = 1:Ne
    x1_A(l,:)=a1_A(l)*s1_A; %x1_A=a1_A*s1_A tin hieu A tai mang
end
x1=x1_A;
x1=awgn(x1,SNR,'measured'); %them nhieu Gaussian trang vao tin hieu
%% (4) Ap dung thuat toan LMS
y = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
e = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
w = zeros(1, Ne); %tao ma tran 0 kich thuoc 1xn1
abswn = zeros(Ne, n1); %tao ma tran 0 kich thuoc Nexn1
wn = zeros(Ne, n1); %tao ma tran 0 kich thuoc Nexn1
for m=1:n1
    y(m)=w*x1(:,m); %dau ra bo loc y
    e(m)=s1_A(m)'-y(m); %sai so/ loi
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
%% (5) Ap dung thuat toan RLS
forgettingfactor=1;
P = (0.5)^(-1)*eye(Ne);
g  = 0;
y1 = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
e1 = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
w1 = zeros(Ne, 1); %tao ma tran 0 kich thuoc Nexn1
P_n = zeros(Ne,Ne);
absw1n = zeros(Ne, n1); %tao ma tran 0 kich thuoc Nexn1
for m = 1:n1
    y1(m)=x1(:,m)'*w1; %dau ra bo loc y
    e1(m)=s1_A(m)-y1(m); %sai so/ loi     
    g=P*x1(:,m)/(forgettingfactor+x1(:,m)'*P*x1(:,m)); %cap nhat gain
    P_n=(P-g*x1(:,m)'*P)/forgettingfactor; %cap nhat ma tran P
    w1=w1+e1(m)*g;% cap nhat trong so w
   	P = P_n;
    absw1n(:,m)=abs(w1); %gia tri abs trong so w
end
a2=zeros(Ne,1);
RLS=zeros(1,length(phi));
for the=1:length(phi) %quet phi tu -180:180
    for l=1:Ne
        a2(l,1)=exp(1j*k*R*cos(phi(the)*pi/180-2*pi*((l-1)/Ne)));
        %tao ma tran vector lai cua mang
    end
    RLS(the)=abs(w1'*a2); %cong thuc pho dap ung thuat toan RLS
end

figure (7); 
plot(phi,20*log10(LMS/max(LMS))); hold on
plot(phi,20*log10(RLS/max(RLS))); 
grid on;

ee=abs(e);
ee1=abs(e1); 
figure(11);plot(t1,ee); hold on; plot(t1,ee1); grid on;
