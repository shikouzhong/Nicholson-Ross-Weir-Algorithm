[f,s11r,s11i,s21r,s21i,s12r,s12i,s22r,s22i] = importfile('11-4.74.txt',7, 227);
% [f,s11r,s11i,s21r,s21i,s12r,s12i,s22r,s22i] = importfile('22-2.txt',7, 227);

%% 输入样品以及夹具长度
holder = 6.950/1000;
position = 0.000/1000;
sample = 4.74/1000;
% sample = 2.00/1000;

%% 定义常数
c = 3.0e8;
fc = 0; % 截止频率
kc = 2 * pi * fc / c;
k = 2 * pi * f / c;

%% S参数矩阵以及测量参考面移动
s11 = s11r + j * s11i;
s21 = s21r + j * s21i;
s12 = s12r + j * s12i;
s22 = s22r + j * s22i;

% 自动调整样品位置
fitnessfcn = @(x) rotate(s11, s21, s12, s22, holder, x, sample, k, kc);
[position, fval] = ga(fitnessfcn, 1, [], [], [], [], -0.001, 0.001)

% 返回测量参考面移动后的S参数矩阵
[fitnessfcn, newS11, newS21, newS12, newS22] = rotate(s11, s21, s12, s22, holder, position, sample, k, kc);

%% NRW算法计算

[er1, ei1, ur1, ui1] = NRW(newS11, newS21, f, sample, k, kc);

[er2, ei2, ur2, ui2] = NRW(newS22, newS12, f, sample, k, kc);

er = 0.5*(er1+er2);
ei = 0.5*(ei1+ei2);
ur = 0.5*(ur1+ur2);
ui = 0.5*(ui1+ui2);

%% 画图比较
figure;
hold on;
% plot(f,ur_m);
plot(f,ur);
plot(f,ur1);plot(f,ur2);
% plot(f,ui_m);
plot(f,ui);
plot(f,ui1);plot(f,ui2);
% legend('ur_m','ur','ur1','ur2','ui_m','ui','ui1','ui2');

figure;
hold on;
% plot(f,er_m);
plot(f,er);
plot(f,er1);plot(f,er2);
% plot(f,ei_m);
plot(f,ei);
plot(f,ei1);plot(f,ei2);
% legend('er_m','er','er1','er2','ei_m','ei','ei1','ei2');

%% 比较群时延
% phase11 = unwrap(atan2(imag(newS11), real(newS11)));
% phase21 = unwrap(atan2(imag(newS21), real(newS21)));
% phase12 = unwrap(atan2(imag(newS12), real(newS12)));
% phase22 = unwrap(atan2(imag(newS22), real(newS22)));

% figure;
% subplot(2, 2, 1)
% plot(f,phase11)
% subplot(2, 2, 2)
% plot(f,phase21)
% subplot(2, 2, 3)
% plot(f,phase12)
% subplot(2, 2, 4)
% plot(f,phase22)

% 群时延测量值
% a = 11; % a大于1的奇数，代表孔径
% tau_m = zeros(length(f) - a + 1, 1);
% newF = zeros(length(tau_m), 1);
% for k = 1 : length(tau_m)
%     tau_m(k) = -(phase21(k + a -1) - phase21(k)) / (f(k + a -1) - f(k)) / (2 * pi);
%     newF(k) = f((2 * k + a - 1) / 2);
% end
% 群时延计算值
% Q = (mu .* epsilon .* f.^2 / c^2 - fc^2 / c^2).^(1 / 2);
% tau_c = zeros(length(f) - a + 1, 1);
% for k = 1 : length(tau_m)
%     tau_c(k) = sample * (Q(k + a -1) - Q(k)) / (f(k + a -1) - f(k));
% end
% 微分曲线找突变
% delta_tau = tau_m - tau_c;
% s_delta = smooth(delta_tau, 5);
% V = diff(delta_tau) ./ diff(newF);
% B = diff(s_delta) ./ diff(newF);

% figure;
% hold on;
% plot(newF, tau_m);
% plot(newF, tau_c);
% plot(newF, tau_c2);
% plot(newF, tau_m - tau_c);
% figure;
% plot(newF(1 : length(V)), V);
% plot(newF(1 : length(B)), B);
