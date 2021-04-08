function colorScatterQQ(X,Y,NBINS)
%figure
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

hold on

xSort = sort(X);
ySort = sort(Y);
plot(xSort(ceil(numel(xSort)*.5)),ySort(ceil(numel(ySort)*.5)),'r.','MarkerSize',15)
text(xSort(ceil(numel(xSort)*.5)),ySort(ceil(numel(ySort)*.5)),'50%')
plot(xSort(ceil(numel(xSort)*.9)),ySort(ceil(numel(ySort)*.9)),'r.','MarkerSize',15)
text(xSort(ceil(numel(xSort)*.9)),ySort(ceil(numel(ySort)*.9)),'90%')
plot(xSort(ceil(numel(xSort)*.99)),ySort(ceil(numel(ySort)*.99)),'r.','MarkerSize',15)
text(xSort(ceil(numel(xSort)*.99)),ySort(ceil(numel(ySort)*.99)),'99%')
plot(xSort(ceil(numel(xSort)*.999)),ySort(ceil(numel(ySort)*.999)),'r.','MarkerSize',15)
text(xSort(ceil(numel(xSort)*.999)),ySort(ceil(numel(ySort)*.999)),'99.9%')


shg