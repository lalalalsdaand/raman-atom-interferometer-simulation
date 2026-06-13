function dydt=free_fall(~, y, gravity)
dydt=[y(4:6); gravity(1:3);zeros(size(y,1)-6,1)];
%原子处于重力场中的函数