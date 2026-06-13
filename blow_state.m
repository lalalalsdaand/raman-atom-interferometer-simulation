function [state] = blow_state(state,F)
%UNTITLED3 此处显示有关此函数的摘要
%

if F==2
    state.rho(1:5,:,:)=zeros(5,8,state.N);
    state.rho(:,1:5,:)=zeros(8,5,state.N);
elseif F==1
    state.rho(6:8,:,:)=zeros(3,8,state.N);
    state.rho(:,6:8,:)=zeros(8,3,state.N);
end
end

