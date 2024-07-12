global func_time func_current;
global func_time_half func_current_half func_conti;
global func_x0;

background_lsq=current;
faradic_lsq=current;
temp_time=time(nidx(1):nidx(2)-1)-time(nidx(1));
background_lsq_half=current;
faradic_lsq_half=current;




for i=1:length(nidx)-1
    temp_i=current(nidx(i):nidx(i+1)-1);
    
    func_time=temp_time;
    func_current=temp_i;
    x0=abs(ft_i(2,i));
    p0=angle(ft_i(2,i));
    o0=mean(temp_i);


    fun = @(A)A(1)*x0*cos(2*pi*f(2)*temp_time+A(2))+A(3)-temp_i;

    func_x0=x0;

    options = optimset('Display','iter');
   
    x=lsqnonlin(fun,[0.2, p0, o0],[0.001,-2*pi,-0.2],[1, 2*pi, 0.2],[],[],[],[],@nonlcon);
    if rmse_value>0.005
        background_lsq(nidx(i):nidx(i+1)-1)=x(1)*x0*cos(2*pi*f(2)*temp_time+x(2))+x(3);
        faradic_lsq(nidx(i):nidx(i+1)-1)=current((nidx(i):nidx(i+1)-1))-background_lsq(nidx(i):nidx(i+1)-1);
    else
        background_lsq(nidx(i):nidx(i+1)-1)=background(nidx(i):nidx(i+1)-1);
        faradic_lsq(nidx(i):nidx(i+1)-1)=current((nidx(i):nidx(i+1)-1))-background(nidx(i):nidx(i+1)-1);
    end

    temp_i_1st=temp_i(1:40);
    temp_i_2nd=temp_i(41:end);
    temp_time_1st=temp_time(1:40);
    func_time_half=temp_time_1st;

    fun_1st = @(A)A(1)*cos(2*pi*f(2)*temp_time_1st+A(2))+A(3)-temp_i_1st;
    func_current_half=temp_i_1st;
   

    func_current_half=temp_i_2nd;
    
    fun_2nd = @(A)A(1)*cos(2*pi*f(2)*temp_time_1st+A(2))+A(3)-temp_i_2nd;
end

%%
figure(100);
clf;
plot(time, current,'k.-',time, background_lsq,'b.-');

figure(101);
clf;
kk=cv_idx;
show_n=1;
plot(potential(nidx(kk):nidx(kk+show_n)-1), current(nidx(kk):nidx(kk+show_n)-1),'k.-');

axis([-0.2 0.6 -0.35 0.35]);
xticklabels({});
yticklabels({});

figure(101);
clf;
kk=cv_idx;
show_n=1;
plot(abs(ft_i(:,kk)),'k.-');
axis([0 20 0 11]);
xticklabels({});
yticklabels({});



figure(102);
clf;
plot(time, current,'k.-',time, background_lsq_half,'b.-');

figure(103);
clf;
kk=cv_idx;
show_n=1;
plot(potential(nidx(kk):nidx(kk+show_n)-1), current(nidx(kk):nidx(kk+show_n)-1),'k.-');
hold on;
plot(potential(nidx(kk):nidx(kk+show_n)-1), background_lsq_half(nidx(kk):nidx(kk+show_n)-1),'b.-');
hold on;
plot(potential(nidx(kk):nidx(kk+show_n)-1), faradic_lsq_half(nidx(kk):nidx(kk+show_n)-1),'r.-');

%%
function [c,ceq] = nonlcon(x)
    global func_time func_current func_x0;
    c = -abs(func_current)+abs(x(1)*func_x0*cos(2*pi*125*func_time+x(2))+x(3));
    %ceq=temp_i_1st+temp_time_1st-x(1);
    ceq=[];
end

function [c,ceq] = nonlcon_1st(x)
    global func_time_half func_current_half;
    c = -func_current_half+(x(1)*cos(2*pi*125*func_time_half+x(2))+x(3));
    %ceq=temp_i_1st+temp_time_1st-x(1);
    ceq=[];
end

function [c,ceq] = nonlcon_2nd(x)
    global func_time_half func_current_half func_conti;
    c = +func_current_half-(x(1)*cos(2*pi*125*func_time_half+x(2))+x(3));
    %ceq=func_conti-(x(1)*cos(2*pi*125*func_time_half(1)+x(2))+x(3));
    ceq=[];
end
