function [fitnessfcn, newS11, newS21, newS12, newS22] = rotate(s11, s21, s12, s22, holder, position, sample, k, kc)
% 变换S参数使得校准参考面移动至MUT两端
%   变换后验证两侧反射相位是否相等

theta_1 = (k.^2 - kc.^2).^(1/2) * position;
theta_2 = (k.^2 - kc.^2).^(1/2) * (holder - position - sample);

newS11 = s11 .* exp(j * 2 * theta_1);
newS21 = s21 .* exp(j * (theta_1 + theta_2));
newS12 = s12 .* exp(j * (theta_1 + theta_2));
newS22 = s22 .* exp(j * 2 * theta_2);

phase_1 = unwrap(atan2(imag(newS11), real(newS11)));
phase_2 = unwrap(atan2(imag(newS22), real(newS22)));
delta_phase = phase_1 - phase_2;
fitnessfcn = sum(delta_phase.^2);

end