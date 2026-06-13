function [state2, tt, yy]=state_eval_3(state,func,tstart, tend, n)
x=state.x;
x2=state.x;
v=state.v;
v2=state.v;
rho=state.rho;
rho2=state.rho;
tt=linspace(tstart,tend,n);
yy=zeros(size(tt,2),3142,state.N);
for ii=1:state.N
%for ii=1:state.N
    [t,y]=ode45(func, [tstart, tend],vertcat(x(:,ii),v(:,ii),reshape(rho(:,:,ii),56*56,1)));% 用ode45函数求解微分方程

    yy(:,:,ii)=interp1(t,y(:,1:3142),tt);%插值
    y=y(end,:);
    x2(:,ii)=y(1:3);
    v2(:,ii)=y(4:6);
    rho2(:,:,ii)=reshape(y(7:end),56,56);
   % rho2(:,:,ii)=extract_rho8_from_blocks(rho40_end);
%     for kk = 1:length(tt)
%         yk = interp1(t,y,tt(kk));
%         rho40_hist(kk,:,:,ii) = reshape(yk(7:1606),40,40);
%     end
end
state2=state;
state2.x=x2;
state2.v=v2;
state2.rho=rho2;


