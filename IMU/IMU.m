%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%姿态解算 Mayhony 'DCM filter' 
%四元数更新为一阶龙库法
%实验效果：在小范围内，能够很好的推算姿态角，大范围内，连续运动，误差最大达到10度
%          原因是该算法并没做全姿态，没有地磁计的数据
%建议可用PID调试，做简单的飞行，但要高精度，不采用
%实验数据不包括地磁数据，因此对Yaw不作初始化，设为0，后期假如，需要改动
%对加速度计数据进行了滤波处理
% 修正了一些代码
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
Ts = 20e-3;
% K_Acc=0.0236;
% K_Gyro=0.0033054;
halfT=0.5*Ts;
Kp=50*Ts;
Ki=10*Ts;
exInt = 0;								 
eyInt = 0;
ezInt = 0; 

fname  = 'ttest02';
evalcmd = [ 'load -ascii ' fname '.dat' ];
eval(evalcmd);
evalcmd = [ 'in_data =' fname ';'];
eval(evalcmd);
Rollrate  = (in_data(:,8))/1000;
Pitchrate = (in_data(:,9))/1000;
Yawrate   = (in_data(:,10))/1000;
Xacc  = (in_data(:,4))/100;
Yacc  = (in_data(:,5))/100;
Zacc  = (in_data(:,6))/100;
gx = Rollrate;        %Roll_rate
gy = Pitchrate;        %Pitch_rate
gz = Yawrate;        %Yaw_rate
[b,a]=butter(2,0.08);  %设计一个滤波器，给加速度的值进行滤波
% ax=filter(b,a,Xacc);   %一维数字滤波器
% ay=filter(b,a,Yacc);
% az=filter(b,a,Zacc);
ax=Xacc;
ay=Yacc;
az=Zacc;
roll_ref=in_data(:,11)/100;       
pitch_ref=in_data(:,12)/100;	   
yaw_ref=in_data(:,13)/100;
[m,n] = size(gx);       

t = (0:(m-1))*Ts;
t_fin = (m-1)*Ts;

q0 = 1;
q1 = 0;
q2 = 0;
q3 = 0;


for i=1:m
    
  % 对四元素进行一个初始化  
  q0q0=q0*q0;
  q0q1=q0*q1;
  q0q2=q0*q2;
  q0q3=q0*q3;
  q1q1=q1*q1;
  q1q2=q1*q2;
  q1q3=q1*q3;
  q2q2=q2*q2;
  q2q3=q2*q3;
  q3q3=q3*q3;

  norm=sqrt(ax(i)*ax(i)+ay(i)*ay(i)+az(i)*az(i));%把加计的三维向量转成单位向量
  ax(i)=ax(i)/norm;
  ay(i)=ay(i)/norm;
  az(i)=az(i)/norm;

  vx = 2*(q1q3 - q0q2);%vx\y\z，其实就是当前的欧拉角（即四元数）的机体坐标参照系上，换算出来的重力单位向量												
  vy = 2*(q0q1 + q2q3);
  vz = q0q0 - q1q1 - q2q2 + q3q3 ;

  ex = (ay(i)*vz - az(i)*vy) ;   %exyz就是两个重力向量的叉积                        					
  ey= (az(i)*vx - ax(i)*vz) ;
  ez = (ax(i)*vy - ay(i)*vx) ;

  exInt = (exInt + ex)* Ki;								 
  eyInt = (eyInt + ey)* Ki;
  ezInt = (ezInt + ez)* Ki; 

  gx(i) = gx(i) + Kp*ex + exInt;					   							
  gy(i) = gy(i) + Kp*ey + eyInt;
  gz(i) = gz(i) + Kp*ez + ezInt;	     

  q0 = q0 + (-q1*gx(i)  - q2*gy(i) - q3*gz(i))*halfT;%四元数微分方程
  q1 = q1 + ( q0*gx(i) + q2*gz(i) - q3*gy(i))*halfT;
  q2 = q2 + ( q0*gy(i)  - q1*gz(i) + q3*gx(i))*halfT;
  q3 = q3 + ( q0*gz(i) + q1*gy(i)  - q2*gx(i))*halfT;

  norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);%四元数规范化
  q0 = q0 / norm;
  q1 = q1 / norm;
  q2 = q2 / norm;
  q3 = q3 / norm;

   Roll(i)     = atan2(2 * q2 * q3 + 2 * q0 * q1, -2 * q1 * q1 - 2 *q2* q2 + 1)* 57.3; %roll
   Pitch(i)   = (asin(2 * q0 * q2 - 2 * q1* q3))* 57.3; % pitch
   Yaw(i)    = atan2(2 * q1 * q2 + 2 * q0 * q3 , -2 *q2 * q2 - 2 * q3 * q3 +1)* 57.3; %yaw
end

subplot(3,1,1)
plot(t,Roll,t,roll_ref,'r');
legend('AHRS caculated','参考 Roll')
subplot(3,1,2)
plot(t,Pitch,t,pitch_ref,'r');
legend('AHRS caculated','参考 Pitch')
subplot(3,1,3)
plot(t,Yaw,t,yaw_ref,'r');
legend('AHRS caculated','参考 Yaw')