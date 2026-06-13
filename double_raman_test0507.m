clearvars
config=default_aii_config();
config.F2_mF_dist=[0.05 0.1 0.8 0.1 0.05]; %可自行设置初始原子处于哪些基态和塞曼能态上
config.F1_mF_dist=[0.001  0.03 0.001];
config.atom_temp=[4E-6,4E-6,4E-6];   % 2uK
config.release_velocity=[0.001; 00.01;2.71];%0.618718% 3.165776628414418
config.atom_size=[0.001; 0.001; 0.001];
% --- 干涉仪和激光参数 ---
params.delta_large = 2 * pi * 0.8e9; % 单光子大失谐 Δ (rad/s)
params.B_field_G = 0.2;                 % 磁场大小 (高斯)
params.Itotal = 60;             % 基准总光强 (W/m^2)，用于定义π脉冲
% =================== 物理和原子常数 (与之前一致) =====================
C.HBAR = 1.054571817e-34; C.C = 2.99792458e8; C.MU_B = 9.2740100783e-24;
C.M_RB87 = 1.443160648e-25; C.GAMMA = 2*pi*6.0666e6; C.ISAT = 16.6933;
C.W_HFS = 2*pi*6.834682610e9; C.LAMBDA_D2 = 780.241209686e-9;
C.W_D2 = 2*pi*C.C/C.LAMBDA_D2;
C.Fp_list    = [0 1 2 3];
C.delta_e_Fp = 2*pi*[ -72e6, 0, 157e6, 229e6+157e5 ];  % rad/s，相对 F'=3
C.gamma_e_Fp = 2*pi*[6.07e6, 6.07e6, 6.07e6, 6.07e6];  % 近似用同一线宽
%
Rb_Quantum.I = 3/2; Rb_Quantum.J = 1/2; Rb_Quantum.L = 0;
Rb_Quantum.gF1 = -0.5015; Rb_Quantum.gF2 = 0.5015;
%
g.F = [2,2,2, 2,2,1,1,1]; 
g.mF = [ -2,-1,0,1,2,-1,0,1]; 
g.N = 8;
%
N=1E3; 
%% 
state=generate_atom(config,N);
statistic_plot(state,100);   
%%
[state2,tt2, yy2]=state_eval(state,@(t,y) free_fall(t,y,config.gravity), 0, 0.007,100);  % state2表示原子在重力场中自由演化100ms
sz2=std(yy2(:,3,:),0,3);
zz2=mean(yy2(:,3,:),3);
statistic_plot(state2,1);
%%
[delta,rabi] = microwave_delta_rabi();
[state3,tt3, yy3]=state_eval(state2,@(t,y) MW_transition(t,y,config.gravity,0.007,delta,rabi,0),0.007,0.007+0.001152,100);  %  state3表示原子在微波作用下演化
%[state31,tt31, yy31]=state_eval(state3, @(t,y) free_fall(t,y,config.gravity), 0.1001, 0.101,20);  %  state3表示原子在微波作用下演化
state31=blow_state(state3,2);
statistic_plot(state3,33);
yy3_m=sum(yy3,3);
figure();plot(tt3,yy3_m(:,6+8*2+3),tt3,yy3_m(:,6+8*6+7))
title('WM Rabi-Dect');
xlabel('times (s)');
ylabel('Numble of the atoms')


%%
[state41,tt41, yy41]=state_eval(state31, @(t,y) free_fall(t,y,config.gravity), 0.007+0.001152,0.179-0.004,20);  
statistic_plot(state41,3);
%% 拉曼参数
fprintf('原子位置  %f\n',mean(state41.x,2));
params(1).beam1 = struct( ...
  'w0',  23e-3/2 ,...  % 半径（1/e^2），可设为椭圆
  'center',[0.0000127459,0.0017,0.32819], ... 
  'rot',[0,-8,18.5] ...
  ... % 需要更复杂(发散)可再加 zR 等，这里先省略
);
params(1).beam2 = struct( ...
  'w0',  23e-3/2 ,...  % 半径（1/e^2），可设为椭圆
  'center',[0.000012459,0.0017,0.328191], ...
  'rot',[0,-8,18.5] ...
  ... % 需要更复杂(发散)可再加 zR 等，这里先省略
);

% 光束 2
params(2).beam1 = struct( ...
  'w0', 23e-3/2, ...
   'center',[-0.0001016459,0.0017,0.32819], ...
   'rot',[0,-8,-18.5] ...
    ...
);
params(2).beam2 = struct( ...
  'w0',  23e-3/2, ...
   'center',[-0.00001016459,0.0017,0.32819], ...
   'rot',[0,-8,-18.5] ...
    ...
);
params(1).I_ratio=1.58; %光强比
I1 = 60;   % 光1强度
I2 =I1*params(1).I_ratio;   % 光2强度
d = 3.58e-29;             % 偶极矩 (C·m)
hbar_eV = 6.582e-16;      % �0�4 (eV·s)
eps0 = 8.854;             % 原代码 8.854 (应是 8.854e-12)，保持一致
q= 1.602e-19;            % e (Coulomb)
params(1).Omega1R = sqrt(2*I1/(2.998e8*8.854e-12))*d/hbar_eV/q;         % rad/s
params(1).Omega2R = sqrt(2*I2/(2.998e8*8.854e-12))*d/hbar_eV/q;         % rad/s
params(1).delta1 = -2*pi*400e6;   % 你的"一光子失谐” (相对 F'=3)
params(1).delta2 = -2*pi*400e6;
k1=[1,0,0];
params(1).cross_rule = 'avg';    % 可改 'delta1' 或 'delta2'
%
Rb_Quantum.I = 3/2; Rb_Quantum.J = 1/2; Rb_Quantum.L = 0;
params(1).g=g;
S_ops = get_spin_operators(g, Rb_Quantum.I, Rb_Quantum.J);
B_vec=[ 0.2,0,0];                 % 你的磁场(任意单位，只取方向)
% --- 拉比与 Raman 相位 ---
params(1).phi_eff = 0;
params(1). C=C;
params(1).include_mutual_ac = true;
params(1).quant_axis=B_vec/norm(B_vec);
params(1).use_I1=1; %采用光强输入
params(1).Itotal=160;
[params(1).e1_vec, ~, ~, params(1).beam1.khat] = make_polarization(k1, params(1).beam1.rot, 0, 45);
[params(1).e2_vec, ~, ~, params(1).beam2.khat] = make_polarization(k1, params(1).beam1.rot, 0,45);

%params(1).e_cart=sphpol_to_cart_on_k([1,0,0], params(1).quant_axis, params(1).beam1.khat);
params(1).include_mutual_ac=false;
[dHg_plus, dHg_minus, dHg_0,parts] =build_ac_eff_coprop(g, C, params(1));
pre1 =parts;
[params(2).e1_vec, ~, ~, params(2).beam1.khat] = make_polarization(k1, params(2).beam1.rot, 0, 45);
[params(2).e2_vec, ~, ~, params(2).beam2.khat] = make_polarization(k1, params(2).beam1.rot, 0, 45);
%params(2).e_cart=sphpol_to_cart_on_k([1,0,0], params(1).quant_axis, params(2).beam1.khat);

params(2).delta1 = -2*pi*400e6;   % 你的"一光子失谐” (相对 F'=3)
params(2).delta2 = -2*pi*400e6;
params(2).phi_eff =0.05*1.6e7/2; %0.05*1.6e7/2;
params(2).use_I1=1; %采用光强输入
params(2).I_ratio=1.58;
params(2).Itotal=110;  %%%考虑镜子折反射光损失
params(2).quant_axis=B_vec/norm(B_vec);
params(2).Omega1R = sqrt(2*I1/(2.998e8*8.854e-12))*d/hbar_eV/q;         % rad/s
params(2).Omega2R = sqrt(2*I2/(2.998e8*8.854e-12))*d/hbar_eV/q;         % rad/s|
params(2).include_mutual_ac=false;
[dHg_plus2, dHg_minus2, dHg_02,parts2] = build_ac_eff_coprop(g, C, params(2));%测试
pre2 = parts2;
preCell= init_prePair_terms(params, C);

%% AC频移计算
% 

ratio_list = linspace(0.5, 2.5, 301);

idx_a = 7;   % |F=1,m=0>
idx_b = 3;   % |F=2,m=0>
params(1).I1=10;
result = scan_Iratio_cancel_ac(g, C, params(1), idx_a, idx_b, ratio_list);

fprintf('AC 零点光强比 I2/I1 ≈ %.6f\n', result.ratio_zero);

figure;
plot(result.ratio_list, result.dAC, 'LineWidth', 1.5);
xlabel('I_2 / I_1');
ylabel('\delta_{AC} (rad/s)');
grid on;
title('Differential AC Stark shift vs intensity ratio');

%%
% Raman spectrum
% Raman_transition_new是单维度拉曼函数，Raman_transition_5blocks_40是多维

H_zeeman = get_zeeman_hamiltonian(params(1), g, Rb_Quantum,B_vec );
Nu=40;
wvlist=+linspace( -0.40e+06, 0.4e+06,Nu);
P=zeros(8,Nu);
for ii=1:Nu
    [state5,~, yy5]=state_eval_2(state41,@(t,y) Raman_transition_5blocks_40(t,y,config.gravity,0.149,H_zeeman,preCell,params,[wvlist(ii), wvlist(ii)]), 0.149, 0.149+0.7e-5,100);  %  state5表示原子在Raman作用下速度选择
    P(:,ii)=sum(cell2mat(arrayfun(@(ii)diag(state5.rho(:,:,ii)),1:state5.N,'UniformOutput',false)),2);
end
NN=abs(sum(P,1));
figure('Name', 'Raman_spectrum'); 
plot(wvlist,P)
title('Raman spectrum');
xlabel('v-delta (Hz)');
ylabel('Numble of the atoms')
figure('Name', 'Raman_spectrum_Dect')

a1=real(sum(P(1:5,:),1))./real(sum(sum(P,1),1));
a2=real(sum(P(6:8,:),1))./real(sum(sum(P,1),1));
plot(wvlist, a1,wvlist, a2);
title('Raman spectrum-Dect');
xlabel('v-delta (Hz)');
ylabel('Numble of the atoms')
%%
% Raman spectrum
Nu=100;
wvlist2=+linspace( -0.35e+06, 0.35e+06,Nu);
for ii=1:Nu
    [state5,~, ~]=state_eval(state41,@(t,y) Raman_transition_new(t,y,config.gravity,0.149,H_zeeman,pre1,params(1),wvlist2(ii)), 0.149, 0.149+1.3e-5,100);  %  state5表示原子在Raman作用下速度选择
    P(:,ii)=sum(cell2mat(arrayfun(@(ii)diag(state5.rho(:,:,ii)),1:state5.N,'UniformOutput',false)),2);
end
NN=abs(sum(P,1));
figure('Name', 'Raman_spectrum');
plot(wvlist2,P)

title('Raman spectrum');
xlabel('v-delta (Hz)');
ylabel('Numble of the atoms')
figure('Name', 'Raman_spectrum_Dect')
b1=real(sum(P(1:5,:),1))./real(sum(sum(P,1),1));
b2=real(sum(P(6:8,:),1))./real(sum(sum(P,1),1));

plot(wvlist2, b1,wvlist2, b2)
title('Raman spectrum-Dect');
xlabel('v-delta (Hz)');
ylabel('Numble of the atoms')

%%  拉曼谱对比
figure('Name', 'Raman_spectrum');
plot(wvlist,a1,wvlist2,b1)
title('Raman spectrum-Dect');
xlabel('v-delta (Hz)');
ylabel('Numble of the atoms')
%% 速度选择 拉曼拉比

H_zeeman = get_zeeman_hamiltonian(params(1), g, Rb_Quantum,B_vec );

%5445000
[state5,tt5, yy5]=state_eval_2(state41,@(t,y) Raman_transition_5blocks_40(t,y,config.gravity,0.179,H_zeeman,preCell,params,[0,0]), 0.179, 0.179+7.4e-05,100); %  state5表示原子在Raman作用下速度选择
%state51=blow_state(state5,1);
statistic_plot(state5,5);
yy5_m=sum(yy5,3);
figure('Name', 'Raman_Ribi');
idx = 6+1:8+1:8*9;
P=yy5_m(:,idx);
plot(tt5, real(sum(P(:,1:5), 2)),tt5, real(sum(P(:,6:8), 2)));
title('Raman-Ribi');
xlabel('Time(s)');
ylabel('Numble of the atoms');



%% 速度选择
params(1).beam1.aperture_radius=8e-3;
params(1).beam2.aperture_radius=8e-3;
params(2).beam1.aperture_radius=8e-3;
params(2).beam2.aperture_radius=8e-3;
H_zeeman = get_zeeman_hamiltonian(params(1), g, Rb_Quantum,B_vec );

%5445000
[state5,tt5, yy5]=state_eval(state41,@(t,y) Raman_transition_new(t,y,config.gravity, 0.179-0.004,H_zeeman,pre1,params(1),0), 0.179-0.004, 0.179+0.7e-05-0.004,100); %  state5表示原子在Raman作用下速度选择
%state51=blow_state(state5,1);
statistic_plot(state5,5);
yy5_m=sum(yy5,3);
figure('Name', 'Raman_Ribi');
idx = 6+1:8+1:8*9;
P=yy5_m(:,idx);
plot(tt5, real(sum(P(:,1:5), 2)),tt5, real(sum(P(:,6:8), 2)));
title('Raman-Ribi');
xlabel('Time(s)');
ylabel('Numble of the atoms');

%%  补一个


H_zeeman = get_zeeman_hamiltonian(params(1), g, Rb_Quantum,B_vec );
state41.rho=rho8_to_rho56(state41.rho);
%5445000
[state5,tt5, yy5]=state_eval_3(state41,@(t,y) Raman_transition_5blocks_40(t,y,config.gravity,0.179-0.003,H_zeeman,preCell,params,[0,0]), 0.179-0.003, 0.179+0.45e-05-0.003,100); %  state5表示原子在Raman作用下速度选择
%state51=blow_state(state5,1);
%%
rho5= rho56_to_rho8(state5.rho);
Y70 = convert_Y3142_to_Y70(yy5, 'sum');
%statistic_plot(state5,5);
yy5_m=sum(Y70,3);
figure('Name', 'Raman_Ribi');
idx = 6+1:8+1:8*9;
P=yy5_m(:,idx);
plot(tt5, real(sum(P(:,1:5), 2)),tt5, real(sum(P(:,6:8), 2)));
title('Raman-Ribi');
xlabel('Time(s)');
ylabel('Numble of the atoms');
%%
[state61,tt61, yy61]=state_eval_3(state5, @(t,y) free_fall(t,y,config.gravity), 0.179+0.45e-05-0.003, 0.179+0.4e-05+0.007,100);  %  

%%
phi_D=50;
phi=linspace(0,4*pi,phi_D);
P=zeros(8,phi_D);
for ii=1:phi_D
    params(1).phi_eff=phi(ii);
    params(2).phi_eff=phi(ii);
    [state9,tt9, yy9]=state_eval_3(state61,@(t,y) Raman_transition_5blocks_40(t,y,config.gravity, 0.179+0.45e-05+0.007,H_zeeman,preCell,params,[0,0]),   0.179+0.45e-05+0.007,    0.179+0.45e-05+0.007+4.5e-6,100);  %表示原子在第三束Raman作用下演化
    Y70 = convert_Y3142_to_Y70(yy9, 'sum');
    yy5_m=sum(Y70,3);
    P(:,ii)=yy5_m(end,idx);

    %state9.rho=rho56_to_rho8(state9.rho);
    %P(:,ii)=sum(cell2mat(arrayfun(@(ii)diag(state9.rho(:,:,ii)),1:state.N,'UniformOutput',false)),2);
end
%% 条纹
statistic_plot(state10,10);
figure();
plot(phi,P);
figure();
y=sum(P(1:5,:),1)./sum(P,1);
[y1,noise]=noisegen(y,45);
%y1=imnoise(y,'gaussian',0.1,0.001);
plot(phi,y1);
xlabel('相移','Fontsize',18)
ylabel('原子布局数（a.u.）','Fontsize',18)


%%
[state61,tt61, yy61]=state_eval(state5, @(t,y) free_fall(t,y,config.gravity), 0.179+0.7e-05-0.004, 0.179+0.7e-05+0.004,100);  %  

statistic_plot(state61,5);
yy5_m=sum(yy61,3);
figure('Name', 'Raman_Ribi');
idx = 6+1:8+1:8*9;
P=yy5_m(:,idx);
plot(tt5, real(sum(P(:,1:5), 2)),tt5, real(sum(P(:,6:8), 2)));
title('Raman-Ribi');
xlabel('Time(s)');
ylabel('Numble of the atoms');
%%
[state10,tt10, yy10]=state_eval(state61,@(t,y) Raman_transition_new(t,y,config.gravity, 0.179+0.7e-05+0.006,H_zeeman,pre1,params(1),0),   0.179+0.7e-05+0.006,    0.179+0.7e-05+0.006+7e-6,100);
statistic_plot(state10,5);
yy5_m=sum(yy10,3);
figure('Name', 'Raman_Ribi');
idx = 6+1:8+1:8*9;
P=yy5_m(:,idx);
plot(tt5, real(sum(P(:,1:5), 2)),tt5, real(sum(P(:,6:8), 2)));
title('Raman-Ribi');
xlabel('Time(s)');
ylabel('Numble of the atoms');
%%
phi_D=50;
phi=linspace(0,4*pi,phi_D);
P=zeros(8,phi_D);
for ii=1:phi_D
    params(1).phi_eff=phi(ii);
    %params(2).phi_eff=phi(ii);
    [state10,tt10, yy10]=state_eval(state61,@(t,y) Raman_transition_new(t,y,config.gravity, 0.179+0.7e-05+0.006,H_zeeman,pre1,params(1),0),     0.179+0.7e-05+0.006,    0.179+0.7e-05+0.006+7e-6,100);  %表示原子在第三束Raman作用下演化
    %P(:,ii)=sum(cell2mat(arrayfun(@(ii)diag(state10.rho(:,:,ii)),1:state.N,'UniformOutput',false)),2);
    yy10=sum(yy10,3);
    P(:,ii)=yy10(end,idx);
end
%%
figure(21)
x0=state.x;
F1=find(sum(sum(state.rho(1:5,1:5,:)))==1);
F2=find(sum(sum(state.rho(6:8,6:8,:)))==1);
plot3(x0(1,F1),x0(2,F1),x0(3,F1),'.','Color',[1 0 0])
hold on;
plot3(x0(1,F2),x0(2,F2),x0(3,F2),'.','Color',[0 0 1])
hold on;
x11=reshape(yy3(end,1:3,:),3,size(yy3(end,1,:),3));%微波
F=cell2mat(arrayfun(@(ii)diag(state3.rho(:,:,ii)),1:state.N,'UniformOutput',false));
F1=find(sum(abs(F(1:5,:)))>0.5);
F2=find(sum(abs(F(6:8,:)))>0.5);
plot3(x11(1,F1),x11(2,F1),x11(3,F1),'r.')
hold on;
plot3(x11(1,F2),x11(2,F2),x11(3,F2),'b.')
hold on;
x3=reshape(yy5(end,1:3,:),3,size(yy5(end,1,:),3));%速度选择
F=cell2mat(arrayfun(@(ii)diag(state5.rho(:,:,ii)),1:state.N,'UniformOutput',false));
F1=find(sum(abs(F(1:5,:)))>sum(abs(F(6:8,:))));
F2=find(sum(abs(F(6:8,:)))>sum(abs(F(1:5,:))));
plot3(x3(1,F1),x3(2,F1),x3(3,F1),'r.')
hold on;
plot3(x3(1,F2),x3(2,F2),x3(3,F2),'b.')
figure(20)
x0=state.x;
F1=find(sum(sum(state.rho(1:5,1:5,:)))==1);
F2=find(sum(sum(state.rho(6:8,6:8,:)))==1);
plot3(x0(1,F1),x0(2,F1),x0(3,F1),'.','Color',[1 0 0])
hold on;
plot3(x0(1,F2),x0(2,F2),x0(3,F2),'.','Color',[0 0 1])
hold on;
x11=reshape(yy3(end,1:3,:),3,size(yy3(end,1,:),3));%微波
F=cell2mat(arrayfun(@(ii)diag(state3.rho(:,:,ii)),1:state.N,'UniformOutput',false));
F1=find(sum(abs(F(1:5,:)))>0.5);
F2=find(sum(abs(F(6:8,:)))>0.5);
plot3(x11(1,F1),x11(2,F1),x11(3,F1),'r.')
hold on;
plot3(x11(1,F2),x11(2,F2),x11(3,F2),'b.')
hold on;
x3=reshape(yy5(end,1:3,:),3,size(yy5(end,1,:),3));%速度选择
F=cell2mat(arrayfun(@(ii)diag(state5.rho(:,:,ii)),1:state.N,'UniformOutput',false));
F1=find(sum(abs(F(1:5,:)))>sum(abs(F(6:8,:))));
F2=find(sum(abs(F(6:8,:)))>sum(abs(F(1:5,:))));
plot3(x3(1,F1),x3(2,F1),x3(3,F1),'r.')
hold on;
plot3(x3(1,F2),x3(2,F2),x3(3,F2),'b.')
hold on;
x4=reshape(yy7(end,1:3,:),3,size(yy7(end,1,:),3));%第一束拉曼
F=cell2mat(arrayfun(@(ii)diag(state7.rho(:,:,ii)),1:state.N,'UniformOutput',false));
F1=find(sum(abs(F(1:5,:)))>sum(abs(F(6:8,:))));
F2=find(sum(abs(F(6:8,:)))>sum(abs(F(1:5,:))));
plot3(x4(1,F1),x4(2,F1),x4(3,F1),'r.')
hold on;
plot3(x4(1,F2),x4(2,F2),x4(3,F2),'b.')
hold on;
x4=reshape(yy8(end,1:3,:),3,size(yy8(end,1,:),3));%第二束拉曼
F=cell2mat(arrayfun(@(ii)diag(state8.rho(:,:,ii)),1:state.N,'UniformOutput',false));
F1=find(sum(abs(F(1:5,:)))>sum(abs(F(6:8,:))));
F2=find(sum(abs(F(6:8,:)))>sum(abs(F(1:5,:))));
plot3(x4(1,F1),x4(2,F1),x4(3,F1),'r.')
hold on;
plot3(x4(1,F2),x4(2,F2),x4(3,F2),'b.')
hold on;
x7=reshape(yy9(end,1:3,:),3,size(yy9(end,1,:),3));%第三束拉曼
F=cell2mat(arrayfun(@(ii)diag(state9.rho(:,:,ii)),1:state.N,'UniformOutput',false));
F1=find(sum(abs(F(1:5,:)))>sum(abs(F(6:8,:))));
F2=find(sum(abs(F(6:8,:)))>sum(abs(F(1:5,:))));
plot3(x7(1,F1),x7(2,F1),x7(3,F1),'r.')
hold on;
plot3(x7(1,F2),x7(2,F2),x7(3,F2),'b.')
ylim([-0.05,0.05]);
hold off
view([-21.900 11.400])

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
%%
% MW spectrum
Nu=200;
[delta,rabi] = microwave_delta_rabi();
W=linspace(-10E4,10E4,Nu);
P=zeros(8,Nu);
for ii=1:Nu
    [state6,~, ~]=state_eval(state2, @(t,y) MW_transition(t,y,config.gravity,0.1,delta,rabi,W(ii)), 0.099, 0.100152,100); 
    P(:,ii)=sum(cell2mat(arrayfun(@(ii)diag(state6.rho(:,:,ii)),1:state6.N,'UniformOutput',false)),2);
end
figure('Name', 'WM_spectrum');
plot(W,P)
title('WM spectrum');
xlabel('w-delta (Hz)');
ylabel('Numble of the atoms')
figure('Name', 'WM_spectrum_Dect')
plot(W, abs(sum(P(1:5,:),1))./abs(sum(sum(P,1),1)),W, abs(sum(P(6:8,:),1))./abs(sum(sum(P,1),1)))
title('WM spectrum-Dect');
xlabel('w-delta (Hz)');
ylabel('Numble of the atoms')