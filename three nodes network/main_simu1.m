%numerical simulation in one and two strong connection cases
clear

N = 8;
out = [];
for i = 1:N
    epsilon = 2^(-(i-1));
    temp = scalecmi(epsilon);
    temp = temp(:);
    out = [out, temp];
end

% figure; plot( 0:(N-1), out(4,:), 'r', 0:(N-1), out(3,:), 'b--',0:(N-1), out(5,:), 'g--', 0:(N-1), out(1,:), 'k--', 'LineWidth',2, 'MarkerSize', 6);
% axis([0,N-1,0,0.2]);
% xlabel('log_2(1/epsilon)','FontSize', 10, 'FontWeight', 'Bold');
% ylabel('index value','FontSize',10,'FontWeight', 'Bold');
% %title('one strong connection', 'FontWeight', 'Bold', 'FontSize', 10)
% title('two strong connections', 'FontWeight', 'Bold', 'FontSize', 10)
% legend('MPMI', 'NPA', 'PMI', 'CMI', 'Location' ,'SouthWest');
% legend('boxoff');
% set(gca, 'FontSize', 10, 'FontWeight', 'Bold')

 figure;
plot(0:(N-1), log2(out(4,:)), 'r', 0:(N-1), log2(out(3,:)), 'b--', 0:(N-1), log2(out(5,:)), 'g--', 0:(N-1), log2(out(1,:)), 'k--', 'LineWidth',2, 'MarkerSize', 6);
xlabel('log_2(1/epsilon)','FontSize', 10, 'FontWeight', 'Bold');
ylabel('log_2(index value)','FontSize',10,'FontWeight', 'Bold');
%title('one strong connection', 'FontWeight', 'Bold', 'FontSize', 10)
title('two strong connections', 'FontWeight', 'Bold', 'FontSize', 10)
legend('MPMI', 'NPA', 'PMI', 'CMI', 'Location' ,'SouthWest');
legend('boxoff');
set(gca, 'FontSize', 10, 'FontWeight', 'Bold')
