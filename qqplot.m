function qqplot(altHs,mdHs)
altMax= max(altHs);
mdMax = max(mdHs);
hsMax = max([altMax mdMax]);

figure
altS = sort(altHs);
mdS  = sort(mdHs);
plot(0:20,0:20,'k--')
hold on 
plot(altS,mdS,'o')
title('q-q')
xlabel('Altimeter Wave Height [m]')
ylabel('Model Wave Height [m]')
grid on
axis([0 hsMax 0 hsMax])
axis('square')
fontsize(20,20,20,20)


