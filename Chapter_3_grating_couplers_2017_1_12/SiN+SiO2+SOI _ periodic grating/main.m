clear 
clc
%% 
%% 列举出基因组序列！

parameters = [];

    for Period = 1.25:0.05:1.5
        for h_sio = 1.2:0.05:2.0
            parameters = [parameters; Period h_sio];
        end
    end 

%% 第一步，计算初始的基因效率。
efficiency = [] ;
for  num = 1:size(parameters,1)
     
     model = mphopen('erwei_grating_coupler_2017年1月12日 氮化硅+二氧化硅+SOI - 反向.mph');
     model.param.set('Period', [num2str(parameters(num,1)) '[um]']);   %代入周期
     model.param.set('h_sio', [num2str(parameters(num,2)) '[um]']);   %代入齿高度
     model.study('std1').run;    
     
% 入射的光功率（用于后续的效率计算）
     in = mphglobal(model,'Powerin'); 
% 反射率用 S11 参数、通过率用S31参数显示。
%      R = mphglobal(model,'abs(ewfd.S11)^2'); %相信反射率用S11参数更准确，其他反射的功率
%      T = mphglobal(model,'abs(ewfd.S31)^2'); %同样，用S31参数来代表，还留在波导中光的功率
% 向上辐射的功率
     up = mphglobal(model,'Powerup');
% 获得     
     Efield = mphinterp(model,'ewfd.Ez','dataset','cln1');      
     
     angle = 8;
     overlapp = [];
     period = parameters(num,1);
     for ang = 1:size(angle,2)
      align = - period *((1:100)/5);    
      %align = align(Position);
     for i = 1: size(align,2)
         z = linspace(-10,30,size(Efield,2));
         for j = 1: size(z,2);
             Efib(j) = exp(-(((align(i)+z(j))/5.2).^2))*exp(1i*1.46*2*pi/1.55*sind(angle(ang))*z(j));  % 这个
         end 
        overlap(i) =   (abs(sum (conj(Efield).*Efib)))^2 ./ (sum(abs(Efield).^2)) / sum(abs(Efib).^2) ;        
     end
       overlapp = [overlapp ; overlap];
end   
     
     
      efficiency = [efficiency; max(overlap)*up];   
      fprintf('初始效率 %d\n',max(overlap)*up)
     
 
     
end