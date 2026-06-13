function [state2, tt, yy]=state_eval(state,func,tstart, tend, n)
x=state.x; 
x2=state.x; 
v=state.v; 
v2=state.v; 
rho=state.rho; 
rho2=state.rho;%8*8矩阵 
tt=linspace(tstart,tend,n); 
yy=zeros(size(tt,2),70,state.N);%70*1 

parfor ii=1:state.N %for ii=1:state.N 
    [t,y]=ode45(func, [tstart, tend],vertcat(x(:,ii),v(:,ii),reshape(rho(:,:,ii),64,1)));% 用ode45函数求解微分方程 
    yy(:,:,ii)=interp1(t,y(:,1:70),tt);%插值 
    y=y(end,:); x2(:,ii)=y(1:3); v2(:,ii)=y(4:6); rho2(:,:,ii)=reshape(y(end,7:70),8,8); 
end 
state2=state; 
state2.x=x2; 
state2.v=v2; 
state2.rho=rho2;
end