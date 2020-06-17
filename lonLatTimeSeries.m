function lonLatTimeSeries(lon,lat,time,altHs,mdHs)

figure
subplot(1,3,1)
plot(time,altHs,'.')
hold on
plot(time,mdHs,'.')
grid on
datetick
legend('altimeter','model')
ylabel('Wave Height [m]');


subplot(1,3,2)
plot(lon,altHs,'.')
hold on
plot(lon,mdHs,'.')
grid on
xlabel('Longitude [\circW]');

subplot(1,3,3)
plot(lat,altHs,'.')
hold on
plot(lat,mdHs,'.')
grid on
xlabel('Lattitude [\circN]');

fontsize(20,20,20,20)
packfig(1,3)
