% plots for Dvory, Yang, & Dunham

% load pseudo-3D model results
% interpolated to 0.5 yr time intervals for plotting

load params0123.mat

figure(1),clf
subplot(3,2,1)
plot(x,p2,'b'),axis([-15 15 0 12])
xlabel('x (km)'),ylabel('pressure change (MPa)')
text(-14,10.25,'k = 5\times10^{-14} m^2')
text(-14,8,'w = 45 m')
hold on,plot([-8 -8],[0 2],'k',[-1 -1],[0 2],'k',[8 8],[0 2],'k'),hold off
subplot(3,2,2)
plot(x,D2,'r'),axis([-15 15 0 0.35])
xlabel('x (km)'),ylabel('slip (m)')
text(-14,0.299,'k = 5\times10^{-14} m^2')
text(-14,0.233,'w = 45 m')
hold on,plot([-8 -8],[0 0.058],'k',[-1 -1],[0 0.058],'k',[8 8],[0 0.058],'k'),hold off

%subplot(3,2,3)
%plot(x,p0,'b'),axis([-15 15 0 12])
%xlabel('x (km)'),ylabel('pressure change (MPa)')
%text(-14,10.25,'k = 10^{-13} m^2')
%text(-14,8,'w = 30 m')
%hold on,plot([-8 -8],[0 2],'k',[-1 -1],[0 2],'k',[8 8],[0 2],'k'),hold off
%subplot(3,2,4)
%plot(x,D0,'r'),axis([-15 15 0 0.35])
%xlabel('x (km)'),ylabel('slip (m)')
%text(-14,0.299,'k = 10^{-13} m^2')
%text(-14,0.233,'w = 30 m')
%hold on,plot([-8 -8],[0 0.058],'k',[-1 -1],[0 0.058],'k',[8 8],[0 0.058],'k'),hold off

subplot(3,2,3)
plot(x,p3(:,1:7),'b',x,p3(:,8:end),'k'),axis([-15 15 0 12])
xlabel('x (km)'),ylabel('pressure change (MPa)')
text(-14,10.25,'k = 10^{-13} m^2')
text(-14,8,'w = 30 m')
hold on,plot([-8 -8],[0 2],'k',[-1 -1],[0 2],'k',[8 8],[0 2],'k'),hold off
subplot(3,2,4)
plot(x,D3(:,1:7),'r',x,D3(:,8:end),'k'),axis([-15 15 0 0.35])
xlabel('x (km)'),ylabel('slip (m)')
text(-14,0.299,'k = 10^{-13} m^2')
text(-14,0.233,'w = 30 m')
hold on,plot([-8 -8],[0 0.058],'k',[-1 -1],[0 0.058],'k',[8 8],[0 0.058],'k'),hold off

subplot(3,2,5)
plot(x,p1,'b'),axis([-15 15 0 12])
xlabel('x (km)'),ylabel('pressure change (MPa)')
text(-14,10.25,'k = 5\times10^{-13} m^2')
text(-14,8,'w = 45 m')
hold on,plot([-8 -8],[0 2],'k',[-1 -1],[0 2],'k',[8 8],[0 2],'k'),hold off
subplot(3,2,6)
plot(x,D1,'r'),axis([-15 15 0 0.35])
xlabel('x (km)'),ylabel('slip (m)')
text(-14,0.299,'k = 5\times10^{-13} m^2')
text(-14,0.233,'w = 45 m')
hold on,plot([-8 -8],[0 0.058],'k',[-1 -1],[0 0.058],'k',[8 8],[0 0.058],'k'),hold off

%print -depsc2 fig4
