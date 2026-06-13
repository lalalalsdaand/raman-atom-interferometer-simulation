function [config]=default_aii_config()
config.MOT_center=[0; 0; 0];  %设MOT中心为坐标原点
config.release_velocity=[0.618718; 0; 3.165776628414418];  %原子初始速度，我这里设置上抛速度Vz=4.4m/s
config.gravity=[0;0;  -9.793324];
config.atom_size=[0.002; 0.002; 0.002]; %原子初始直径大小
config.atom_temp=2E-6;    %2uK 原子初始温度

config.F2_mF_dist=[1 1 1 1 1]/5*0.98;  %设原子最初处于F=1，和F=2态的各个塞曼能态mF上。（F=1，mF=-1,0,+1,  F=2,mF=-2,-1,0,+1,+2)
config.F1_mF_dist=[1 1 1 ]/3*0.02;

config.B_field_time=0.03;   %后面这些参数我没用到，没仔细看
config.B_field_duration=0.14;

config.MW1_time=0.005;
config.MW1_duration=4E-4;
config.MW1_Rabi=pi/4.2E-4;

config.Blow2_eff=0.99;
config.Blow2_F1_rate=0.005;

config.Raman1_time=0.007;
config.Raman1_duration=1E-5;
config.Raman1_Rabi=pi/1.03E-5;

config.Blow1_eff=0.99;
config.Blow1_F2_rate=0.005;

config.Raman_center=[0; 0.001; 0];
k0=1.61E7;
config.Raman_k=[sin(1E-5) 0 cos(1E-5)] * k0;

config.Raman1_time=0.007;
config.Raman1_duration=1E-5;
config.Raman1_Rabi=pi/1.03E-5;

config.Raman2_time=0.02;
config.Raman2_duration=5E-6;
config.Raman2_Rabi=pi/1.03E-5;

config.Raman3_time=0.075;
config.Raman3_duration=1E-5;
config.Raman3_Rabi=pi/1.03E-5;

config.Raman4_time=0.130;
config.Raman4_duration=5E-6;
config.Raman4_Rabi=pi/1.03E-5;


