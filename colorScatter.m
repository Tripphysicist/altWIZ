function colorScatter(X,Y,NBINS)
figure
histogram2(X,Y,NBINS,'DisplayStyle','tile')
colormap('lansey')
hold on
plot(0:20,0:20,'k--')
%plot(X,Y,'k.')
axis('square')
maxX = max(X); maxY = max(Y);
maxVal = max(maxX,maxY);
minX = min(X); minY = min(Y);
minVal = max(minX,minY);
axis([minVal maxVal minVal maxVal])
fontsize(20,20,20,20)
xlabel('Altimeter Wave Height [m]')
ylabel('Model Wave Height [m]')
colorbar
grid on
shg