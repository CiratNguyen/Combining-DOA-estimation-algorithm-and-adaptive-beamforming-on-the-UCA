close all
clear
clc

%% (1) Thiết lập mô hình tín hiệu, mô hình mảng UCA 2-D và thông số mô phỏng
Fs = 48000; Ts=1/Fs; % Fsample: tần số (Hz), Tsample: chu kỳ lấy mẫu (s)
t0=Ts:Ts:0.1; % Thời gian lấy mẫu (s)
f1=1325; s0_A=1*sin(2*pi*f1*t0).'; % Signal_A (mV)
f2=1350; s0_B=1*sin(2*pi*f2*t0).'; % Signal_B (mV)
f3=1375; s0_C=1*sin(2*pi*f3*t0).'; % Signal_C (mV)
n0=length(t0); % Snapshots - NoS (Number of Samples) - Độ dài tín hiệu
SNR=30; % Tỷ số tín hiệu trên nhiễu
angles_A=[-120 10]*pi/180; % Góc tới Signal_A [azimuth elevation]
angles_B=[0 40]*pi/180; % Góc tới Signal_B
angles_C=[70 50]*pi/180; % Góc tới Signal_C
c = 343; % Tốc độ lan truyền âm thanh trong không khí
lambda_max=c/min([f1 f2 f3]); % Lambda max
Ne=6; % Số lượng phần tử mảng
D=3; % Số lượng tín hiệu tới
R=0.5*Ne*lambda_max/(2*pi); % Bán kính mảng
k=2*pi/lambda_max; % Hệ số hoạt động của mảng
k1=2*pi/(c/f1); % Hệ số góc mỗi tín hiệu
k2=2*pi/(c/f2); % Hệ số góc mỗi tín hiệu
k3=2*pi/(c/f3); % Hệ số góc mỗi tín hiệu
phi=-180:0.05:180; % Góc quét phi để tạo ma trận lái mảng
theta=0:0.05:90; % Góc quét theta để tạo ma trận lái mảng
a0_A=zeros(Ne,1); % Tạo ma trận 0 kích thước Nex1
a0_B=zeros(Ne,1); % Tạo ma trận 0 kích thước Nex1
a0_C=zeros(Ne,1); % Tạo ma trận 0 kích thước Nex1
for l = 1:Ne
    a0_A(l)=exp(1j*k1*R*sin(angles_A(2))*cos(angles_A(1)-2*pi*((l-1)/Ne)));
    % Vector lái Signal_A
    a0_B(l)=exp(1j*k2*R*sin(angles_B(2))*cos(angles_B(1)-2*pi*((l-1)/Ne)));
    % Vector lái Signal_B
    a0_C(l)=exp(1j*k3*R*sin(angles_C(2))*cos(angles_C(1)-2*pi*((l-1)/Ne)));
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
Pmusic=zeros(length(theta),length(phi));
% Tạo ma trận 0 kích thước length(theta),length(phi)
for tt=1:length(theta) %quet theta tu 0:90
    for pp=1:length(phi) %quet phi tu -180:180
        for l=1:Ne
            a0(l,1)=exp(1j*k*R*sin(theta(tt)*pi/180)*cos(phi(pp)*pi/180-2*pi*((l-1)/Ne)));
            %tao ma tran vector lai cua mang
        end
        Pmusic(tt,pp)=abs(1/(a0'*En*(En')*a0)); %cong thuc thuat toan MUSIC
    end
end
%Bieu dien thuat toan MUSIC
figure(1); title('Thuật toán MUSIC') % Biểu đồ phổ MUSIC với giá trị chuẩn hóa
surf(phi,theta,10*log10(Pmusic/max(Pmusic(:))),'EdgeColor', 'none'); %-> dB voi gia tri max = 0
colormap('jet'); colorbar; xlabel('Azimuth plane'); ylabel('Elevation plane'); grid on;
local_maxima1 = imregionalmax(Pmusic);
% Sử dụng hàm imregionalmax để tìm các điểm cục bộ tối đa trong ma trận Pmusic
[row1, col1] = find(local_maxima1); % Lấy tọa độ của các điểm cục bộ tối đa
[max_values1, indices1] = maxk(Pmusic(local_maxima1), D); % Lấy D điểm cục đại cao nhất
max_row1 = sort(row1(indices1));
max_col1 = sort(col1(indices1));
% Hiển thị kết quả
fprintf('Thuật toán MUSIC\n');
for i = 1:length(max_row1)
    angle_1(i, :) = [theta(max_row1(i)), phi(max_col1(i))];
    fprintf('Góc %d:\n [%f, %f]\n', i, phi(max_col1(i)), theta(max_row1(i)));
end
%% (3) Thiet lap thong so dau vao truoc khi ap dung Beamforming
% [biendoPmusic,vitrigoctoi]=findpeaks(Pmusic); %xac dinh goc toi tin hieu
% DOA_est = theta(vitrigoctoi); %qua thuat toan MUSIC
% angles_A1=DOA_est(1)*pi/180; %xac dinh goc toi tin hieu A
% angles_B1=DOA_est(2)*pi/180; %xac dinh goc toi tin hieu B
% angles_C1=DOA_est(3)*pi/180; %xac dinh goc toi tin hieu C
t1=Ts:Ts:5; %thoi gian lay mau (s)
s1_A=1*sin(2*pi*f1*t1).'; %signal_A (mV)
s1_B=1*sin(2*pi*f2*t1).'; %signal_B (mV)
s1_C=1*sin(2*pi*f3*t1).'; %signal_C (mV)
s1=s1_A+s1_B+s1_C; %mixed signal
% [Pxx, ff1] = pwelch(s1, [], [], [], Fs);
n1=length(t1); %bien thoi gian - snapshots
mu = 0.001; %learning rate/ toc do hoc/ he so thich nghi
a1_A=zeros(Ne,1); %tao ma tran 0 kich thuoc Nex1
a1_B=zeros(Ne,1); %tao ma tran 0 kich thuoc Nex1
a1_C=zeros(Ne,1); %tao ma tran 0 kich thuoc Nex1
for l = 1:Ne
    a1_A(l)=exp(1j*k1*R*sin(angles_A(2))*cos(angles_A(1)-2*pi*((l-1)/Ne))); %vector lai signal_A
    a1_B(l)=exp(1j*k2*R*sin(angles_B(2))*cos(angles_B(1)-2*pi*((l-1)/Ne))); %vector lai signal_B
    a1_C(l)=exp(1j*k3*R*sin(angles_C(2))*cos(angles_C(1)-2*pi*((l-1)/Ne))); %vector lai signal_C
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

% luachonphuongphap = input('Lua chon thuat toan LMS (1) hoac RLS (2): ');
% switch luachonphuongphap
% case 1
%% (4) Ap dung thuat toan LMS
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
LMS=zeros(length(theta),length(phi));
for tt=1:length(theta) %quet theta tu -180:180
    for pp=1:length(phi)
        for l=1:Ne
            a1(l,1)=exp(1j*k*R*sin(theta(tt)*pi/180)*cos(phi(pp)*pi/180-2*pi*((l-1)/Ne)));
            %tao ma tran vector lai cua mang
        end
        LMS(tt,pp)=abs(w*a1); %cong thuc pho dap ung thuat toan LMS
    end
end
%Bieu dien thuat toan LMS
figure (2); 
surf(phi,theta,20*log10(LMS/max(LMS(:))),'EdgeColor', 'none'); %-> dB voi gia tri max = 0
grid on;
colormap('jet'); colorbar; 
% figure (3); 
% subplot(311); 
% plot(t1, s1); 
% axis([0.4 0.45 -5 5]);
% ylabel('Biên độ (mV)'); 
% xlabel('Thời gian (s)'); 
% title('TH hỗn hợp'); 
% grid on;
% 
% subplot(312);
% s1_C=awgn(s1_C,SNR,'measured');
% plot(t1,s1_C) ; 
% axis([0.4 0.45 -2 2]);
% ylabel('Biên độ (mV)'); 
% xlabel('Thời gian (s)'); 
% title('TH mong muốn'); 
% grid on;
% 
% subplot(313); 
% plot(t1,y); 
% axis([0.4 0.45 -2 2]);
% ylabel('Biên độ (mV)'); 
% xlabel('Thời gian (s)'); 
% title('TH đầu ra'); 
% grid on;
% 
% [Pyy, ff2] = pwelch(y, [], [], [], Fs); 
% figure(4); 
% plot(ff1, 10*log10(Pxx)); 
% hold on; 
% plot(ff2, 10*log10(Pyy)); 
% axis([900 1700 -85 10]);
% ylabel('Biên độ (dB)'); 
% xlabel('Tần số (Hz)'); 
% grid on; 
% legend('TH hỗn hợp','TH đầu ra');
% 
% figure(5);
% plot(t1,abswn); 
% xlabel('Thời gian (s)'); 
% ylabel('Giá trị trọng số w'); 
% grid on;
% 
% ee=abs(e); 
% figure(6);plot(t1,ee); 
% xlabel('Thời gian (s)'); 
% ylabel('Lỗi'); 
% grid on; 
% axis([0 5 0 2]);

% case 2
% %% (5) Ap dung thuat toan RLS
% forgettingfactor=1;
% P = (0.5)^(-1)*eye(Ne);
% g  = 0;
% y1 = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
% e1 = zeros(1, n1); %tao ma tran 0 kich thuoc 1xn1
% w1 = zeros(Ne, 1); %tao ma tran 0 kich thuoc Nexn1
% P_n = zeros(Ne,Ne);
% absw1n = zeros(Ne, n1); %tao ma tran 0 kich thuoc Nexn1
% for m = 1:n1
%     y1(m)=x1(:,m)'*w1; %dau ra bo loc y
%     e1(m)=s1_C(m)-y1(m); %sai so/ loi     
%     g=P*x1(:,m)/(forgettingfactor+x1(:,m)'*P*x1(:,m)); %cap nhat gain
%     P_n=(P-g*x1(:,m)'*P)/forgettingfactor; %cap nhat ma tran P
%     w1=w1+e1(m)*g;% cap nhat trong so w
%    	P = P_n;
%     absw1n(:,m)=abs(w1); %gia tri abs trong so w
% end
% a2=zeros(Ne,1);
% RLS=zeros(1,length(theta));
% for the=1:length(theta) %quet theta tu -180:180
%     for l=1:Ne
%         a2(l,1)=exp(1j*k*R*cos(theta(the)*pi/180-2*pi*((l-1)/Ne)));
%         %tao ma tran vector lai cua mang
%     end
%     RLS(the)=abs(w1'*a2); %cong thuc pho dap ung thuat toan RLS
% end
% 
% figure (7); 
% plot(theta,20*log10(RLS/max(RLS))); 
% axis([-180 180 -80 0]);
% xlabel('Theta (degree)'); 
% ylabel('Phổ đáp ứng (dB)'); 
% grid on;
% 
% figure(8) 
% subplot(311); 
% plot(t1, s1); 
% axis([0.4 0.45 -5 5]);
% ylabel('Biên độ (mV)'); 
% xlabel('Thời gian (s)'); 
% title('TH hỗn hợp'); 
% grid on;
% 
% subplot(312);
% s1_C=awgn(s1_C,SNR,'measured');
% plot(t1,s1_C) ; 
% axis([0.4 0.45 -2 2]);
% ylabel('Biên độ (mV)'); 
% xlabel('Thời gian (s)'); 
% title('TH mong muốn'); 
% grid on;
% 
% subplot(313); 
% plot(t1,y1); 
% ylabel('Biên độ (mV)'); 
% xlabel('Thời gian (s)'); 
% title('TH đầu ra'); 
% grid on;
% axis([0.4 0.45 -2 2]);
% 
% [Pyyy, fff2] = pwelch(y1, [], [], [], Fs); 
% figure(9); 
% plot(ff1, 10*log10(Pxx)); 
% hold on; 
% plot(fff2, 10*log10(Pyyy)); 
% axis([900 1700 -85 10]);
% ylabel('Biên độ (dB)'); 
% xlabel('Tần số (Hz)'); 
% grid on; 
% legend('TH hỗn hợp','TH đầu ra');
% 
% figure(10);
% plot(t1,absw1n); 
% xlabel('Thời gian (s)'); 
% ylabel('Giá trị trọng số w'); 
% grid on;
% 
% ee1=abs(e1); 
% figure(11);plot(t1,ee1); 
% xlabel('Thời gian (s)'); 
% ylabel('Lỗi'); 
% grid on; 
% axis([0 5 0 2]);
% 
% otherwise
% disp('Chon lai (1) hoac (2)')
% end
% %}
% % figure(7); polarPmusic = polarplot(deg2rad((0:0.01:360)-180),10*log(P)); hold on;
% % APLMS1=APLMS-min(APLMS); polarplot(deg2rad((0:0.01:360)-180),2.5*APLMS1);
% % thetaticks([0 45 90 135 180 225 270 315]); rticks(0); rticklabels('0'); grid on;
% % thetaticklabels({'0°', '45°', '90°', '135°', '±180°', '-135°', '-90°', '-45°'});
% % legend('P_M_U_S_I_C','P_L_M_S'); title('Tương quan giữa phổ đáp ứng P_M_U_S_I_C và P_L_M_S');
% % figure(8);Nee=1:1:6;surf(Nee,t,ww);colormap(cool);


