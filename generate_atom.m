function [state]=generate_atom(config, N)
state.N=N;  %随机生成N个原子（行向量）
state.x=config.MOT_center+config.atom_size.*randn(3,N);  % 表示原子的x,y,z位置坐标，原子初始尺寸服从标准正态分布
state.v=config.release_velocity+sqrt(config.atom_temp*1.38e-23/1.4447e-25)'.*randn(3,N);  %表示原子x,y,z方向的速度，第二项为原子自由膨胀速度，也服从标准正态分布
prob1=horzcat(config.F2_mF_dist,config.F1_mF_dist); 
prob=cumsum(prob1);
r=rand(1,N);
lsi=(r<=prob');
ratio=double(xor([false(1,N); lsi(1:end-1,:)],lsi(1:end,:)));
rho=arrayfun(@(ii) diag(ratio(:,ii)), 1:N, 'UniformOutput', false);
state.rho=cat(3,rho{:});


