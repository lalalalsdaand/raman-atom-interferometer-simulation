function [prob]=statistic_plot(state,fig_num)
figure( fig_num);
subplot(2,4,1)
plot(state.x(1,:),state.x(2,:),'.')    %生成原子y-x,z-x,z-y坐标图
xlabel('x','Fontsize',14)
ylabel('y','Fontsize',14)

daspect([1 1 1])

subplot(2,4,2)
plot(state.x(1 ,:),state.x(3,:),'.')
xlabel('x','Fontsize',14)
ylabel('z','Fontsize',14)
daspect([1 1 1])

subplot(2,4,3)
plot(state.x(2,:),state.x(3,:),'.')
xlabel('y','Fontsize',14)
ylabel('z','Fontsize',14)
daspect([1 1 1])

subplot(2,4,5)
plot(state.x(1,:),state.v(1,:),'.')   %生成原子Vx-x,Vy-y,Vz-z的动量-位置相空间图
xlabel('x','Fontsize',14)
ylabel('vx','Fontsize',14)
pbaspect([1 1 1])

subplot(2,4,6)
plot(state.x(2,:),state.v(2,:),'.')
xlabel('y','Fontsize',14)
ylabel('vy','Fontsize',14)
pbaspect([1 1 1])

subplot(2,4,7)
plot(state.x(3,:),state.v(3,:),'.')
xlabel('z','Fontsize',14)
ylabel('vz','Fontsize',14)
pbaspect([1 1 1])

prob=cell2mat(arrayfun(@(ii)diag(state.rho(:,:,ii)),1:state.N,'UniformOutput',false));
prob=sum(prob,2);
subplot(2,4,4);
bar(prob(1:5));
xlim([0,6]);
xticklabels({'-2','-1','0','1','2'});
xlabel('F=2,mF分布','Fontsize',14)
subplot(2,4,8);
bar(prob(6:8));
xlim([0,4]);
xticklabels({'-1','0','1'});
xlabel('F=1,mF分布','Fontsize',14)