clear 
clc
%% 
%% 列举出基因组序列！
t1 = 0.475;
g1 = 0.475;
t2 = 0.475;
g2 = 0.475;
t3 = 0.475;
g3 = 0.475;
t4 = 0.475;
g4 = 0.475;
t5 = 0.475;
g5 = 0.475;
t6 = 0.475;
g6 = 0.475;
t7 = 0.475;
g7 = 0.475;
t8 = 0.475;
g8 = 0.475;
t9 = 0.475;
g9 = 0.475;
t10 = 0.475;
g10 = 0.475;
t11 = 0.475;
g11 = 0.475;
t12 = 0.475;
g12 = 0.475;
t13 = 0.475;
h_sio = 1.6;


parameters = [  ...
    t1 g1 t2 g2 t3  g3  t4  g4  t5  g5  t6 g6 t7 g7 ...
    t8 g8 t9 g9 t10 g10 t11 g11 t12 g12 t13 ...
    h_sio];
parametersname = {    ...
    't1' 'g1' 't2' 'g2' 't3'  'g3'  't4'  'g4'  't5'  'g5' 't6' 'g6' 't7' 'g7'...
    't8' 'g8' 't9' 'g9' 't10' 'g10' 't11' 'g11' 't12' 'g12' 't13' ...
    'h_sio' }; % 仅供后续对照，避免出现对应不上的情况
%strjoin



%% 第一步，计算初始的基因效率。
     period = 1;
     wl =1.550;
     
     model = mphopen('erwei_grating_coupler_.mph');
     model.param.set('lambda0', [num2str(wl) '[um]']);   %代入波长
     for i = 1 : size(parameters,2)  % 运用循环，代入所有公式
        model.param.set(strjoin(parametersname(i)), [num2str(parameters(i)) '[um]']); 
     end
     model.study('std1').run;    
     
% 入射的光功率（用于后续的效率计算）
     in = mphglobal(model,'Powerin'); 
% 反射率用 S11 参数、通过率用S31参数显示。
     R = mphglobal(model,'abs(ewfd.S11)^2'); %相信反射率用S11参数更准确，其他反射的功率
     T = mphglobal(model,'abs(ewfd.S31)^2'); %同样，用S31参数来代表，还留在波导中光的功率
% 向上辐射的功率
     up = mphglobal(model,'Powerup');
% 获得     
     Efield = mphinterp(model,'ewfd.Ez','dataset','cln1');      
     
     angle = 8;
     overlapp = [];
     for ang = 1:size(angle,2)
      align = - period *((1:100)/5);    
      %align = align(Position);
     for i = 1: size(align,2)
         z = linspace(-10,30,size(Efield,2));
         for j = 1: size(z,2);
             Efib(j) = exp(-(((align(i)+z(j))/5.2).^2))*exp(1i*1.46*2*pi/wl*sind(angle(ang))*z(j));  % 这个
         end 
        overlap(i) =   (abs(sum (conj(Efield).*Efib)))^2 ./ (sum(abs(Efield).^2)) / sum(abs(Efib).^2) ;        
     end
       overlapp = [overlapp ; overlap];
end   
     
     
      efficiency = [max(overlap)*up];   
      fprintf('初始效率 %d\n',efficiency)%显示当前效率
     
     Position = find(overlap==max(overlap));
     
    
%% 接下来进行遗传算法计算。
flag = true;
while (flag)
     flag = false; % 重置flag，用于
     for i = 1 : size(parameters,2)-1;
         temp = parameters;  % 引入临时基因，用于计算
         for temp2 = [ temp(i)-0.005 , temp(i)+0.005 ]; 
             temp(i) = temp2;
             
             
             %第i个参数微调10nm                      
             model = mphopen('erwei_grating_coupler_.mph');
             model.param.set('lambda0', [num2str(wl) '[um]']);   %代入波长
               for j = 1 : size(parameters,2)  % 运用循环，代入所有公式
                   model.param.set(strjoin(parametersname(j)), [num2str(temp(j)) '[um]']);  
               end 
             model.study('std1').run;
             up = mphglobal(model,'Powerup');
             Efield = mphinterp(model,'ewfd.Ez','dataset','cln1');         
             align = - period *((1:100)/5);    
             overlap = [];
             for k = 1: size(align,2)
                 z = linspace(-10,30,size(Efield,2));
                 Efib= [];
                 for p = 1: size(z,2);
                     Efib(p) = exp(-(((align(k)+z(p))/5.2).^2))*exp(1i*1.46*2*pi/wl*sind(8)*z(p));  % 这个
                 end 
                 overlap(k) =   (abs(sum (conj(Efield).*Efib)))^2 ./ (sum(abs(Efield).^2)) / sum(abs(Efib).^2) ;        
             end
             efftemp = max(overlap)*up; %临时基因效率
             
             if efftemp > efficiency  
                 %如果效率增加了
                efficiency = efftemp; %保留效率
                parameters = temp; %保留基因
                fprintf(' %d\n',efficiency)%显示当前效率
                flag = true; %插flag          
             end 
         end
     end
end