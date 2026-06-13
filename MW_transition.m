function dydt=MW_transition(t,y,gravity, tstrat,delta,rabi,w_v)
format long;
% y1= [y(1:8)';horzcat(zeros(1,1),y(9:15)');horzcat(zeros(1,2),y(16:21)');
%     horzcat(zeros(1,3),y(22:26)');horzcat(zeros(1,4),y(27:30)');horzcat(zeros(1,5),y(31:33)');
%     horzcat(zeros(1,6),y(34:35)');horzcat(zeros(1,7),y(36)')];
% y2=conj(y1');
% y2=y2-diag(diag(y2));
% y=y1+y2;
rho=reshape(y(7:end),8,8);
t=t-tstrat;
w=w_v*2*pi;
weff=exp(1j*(delta-w)*t).*rabi ;
H1=[0 0 0 0 0 weff(1,1) weff(1,2) weff(1,3)
    0 0 0 0 0 weff(2,1) weff(2,2) weff(2,3)
    0 0 0 0 0 weff(3,1) weff(3,2) weff(3,3)
    0 0 0 0 0 weff(4,1) weff(4,2) weff(4,3)
    0 0 0 0 0 weff(5,1) weff(5,2) weff(5,3)
    conj(weff(1,1)) conj(weff(2,1)) conj(weff(3,1)) conj(weff(4,1)) conj(weff(5,1)) 0 0 0 
    conj(weff(1,2)) conj(weff(2,2)) conj(weff(3,2)) conj(weff(4,2)) conj(weff(5,2)) 0 0 0
   conj(weff(1,3)) conj(weff(2,3)) conj(weff(3,3)) conj(weff(4,3)) conj(weff(5,3)) 0 0 0]; 

H=(H1)*2*pi*433;%433
dydt=[y(4:6); gravity(1:3); reshape((H*rho-rho*H)/1i,64,1)];
