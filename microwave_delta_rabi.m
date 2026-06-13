function [delta,rabi] = microwave_delta_rabi()
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%磁场
Bc=0.9e-1; %gauss
Ab=1399000;
J=1/2;
I=3/2;
S=1/2;
F1=1;
F2=2;
L=0;
g_j=1+(J*(J+1)+S*(S+1)-L*(L+1))/2/J/(J+1);
g_1=g_j*((F1*(F1+1)-I*(I+1)+J*(J+1))/(F1+1))/2/F1;
g_2=g_j*((F2*(F2+1)-I*(I+1)+J*(J+1))/(F2+1))/2/F2;



VB11=Bc*Ab*g_1*1;
VB1_1=Bc*Ab*g_1*-1;
VB21=Bc*Ab*g_2*1;
VB22=Bc*Ab*g_2*2;
VB2_1=Bc*Ab*g_2*-1;
VB2_2=Bc*Ab*g_2*-2;
%线性极化
rabi1=[0 0 0
     sqrt(3)/4   0   0
     0 1/2  0
      0 0   sqrt(3)/4
      0 0  0 ];
%非线性极化左
rabil=[1/4 0 0
      0  1/12   0
     0 0  1/12
      0 0  0
      0 0  0 ];
%非线性极化右
rabir=[0 0 0
     0   0   0
     1/12 0  0
      0 1/12   0
      0 0  1/4 ];
rabi=rabi1+rabil+rabir;
Bx=2.5e-1;
By=0.03e-1;
VBX11=Bx*Ab*g_1*1;
VBX1_1=Bx*Ab*g_1*-1;
VBX21=Bx*Ab*g_2*1;
VBX22=Bx*Ab*g_2*2;
VBX2_1=Bx*Ab*g_2*-1;
VBX2_2=Bx*Ab*g_2*-2;
VBY11=By*Ab*g_1*1;
VBY1_1=By*Ab*g_1*-1;
VBY21=By*Ab*g_2*1;
VBY22=By*Ab*g_2*2;
VBY2_1=By*Ab*g_2*-1;
VBY2_2=By*Ab*g_2*-2;
delta=[(VB11/2-VB2_2) 0 0
    (VB11+VB2_1)/2   VB2_1/2    (VB1_1+VB2_1)/2
       VB11/2          0       VB1_1/2
       (VB11+VB21)/2  VB21/2   (VB1_1+VB21)/2
        0 0  (VB1_1/2-VB22)].*2*pi;
end

