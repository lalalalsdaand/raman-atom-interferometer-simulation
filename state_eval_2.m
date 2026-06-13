function [state2, tt, yy]=state_eval_2(state,func,tstart, tend, n)
x=state.x;
x2=state.x;
v=state.v;
v2=state.v;
rho=state.rho;
rho2=state.rho;
tt=linspace(tstart,tend,n);
yy=zeros(size(tt,2),70,state.N);
parfor ii=1:state.N
%for ii=1:state.N
    [t,y]=ode45(func, [tstart, tend],init_blocks_state(x(:,ii),v(:,ii),rho(:,:,ii)));% 用ode45函数求解微分方程

    yy(:,:,ii)=interp1(t,convert_y1606_to_y70(y(:,1:3142)),tt);%插值
    y=convert_y1606_to_y70(y(end,:));
    x2(:,ii)=y(1:3);
    v2(:,ii)=y(4:6);
    rho2(:,:,ii)=reshape(y(7:end),8,8);
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


