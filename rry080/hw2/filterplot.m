function output=filterplot(H)

H=H(1,:);
t1=size(H);
f=linspace(-.5,.5,t1(2));

subplot(3,1,1)
plot(f,20*log10(abs(fftshift(fft(H)))))
xlabel('Uniform frequency')
ylabel('Power (dB)')
title('Filter properties')
grid

subplot(3,1,2)
a=unwrap(angle(fftshift(fft(H))));
plot(f,a)
xlabel('Uniform frequency')
ylabel('Angle (rad)')
axis([-.5 .5 min(a) 0])
grid

subplot(3,1,3)
a=grpdelay(H,1);
plot(linspace(-.5,.5,512),fftshift(a)-t1(2)/2)
xlabel('Uniform frequency')
ylabel('Group delay')
axis([-.5 .5 -500 500])
grid

clear t1 f a
