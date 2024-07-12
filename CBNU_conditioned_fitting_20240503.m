global func_time func_current;
global func_time_half func_current_half func_conti;
global func_x0;

background_lsq=current;
faradic_lsq=current;
temp_time=time(nidx(1):nidx(2)-1)-time(nidx(1));
background_lsq_half=current;
faradic_lsq_half=current;


%background_lsq(1:nidx(571)-1)=background(1:nidx(571)-1);
%faradic_lsq(1:nidx(571)-1)=faradic(1:nidx(571)-1);




for i=1:length(nidx)-1
%for i=1000:1010
    temp_i=current(nidx(i):nidx(i+1)-1);
    %temp_i=smooth(temp_time, temp_i, 0.2, 'rloess');
    
    func_time=temp_time;
    func_current=temp_i;
    x0=abs(ft_i(2,i));
    p0=angle(ft_i(2,i));
    o0=mean(temp_i);


    fun = @(A)A(1)*x0*cos(2*pi*f(2)*temp_time+A(2))+A(3)-temp_i;
    %fun = @(A)A(1)*cos(2*pi*f(2)*temp_time+A(2))+A(3);

    func_x0=x0;

    options = optimset('Display','iter');
    %x=lsqnonlin(fun,[0.2, p0, o0]);
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
    %x_1st=lsqnonlin(fun_1st,[x0, p0, o0]);
    func_current_half=temp_i_1st;
    %x_1st=lsqnonlin(fun_1st,[x0, p0, o0],[],[],[],[],[],[],@nonlcon_1st);

    func_current_half=temp_i_2nd;
    %func_conti=x_1st(1)*cos(2*pi*f(2)*temp_time_1st(end)+x_1st(2))+x_1st(3)
    fun_2nd = @(A)A(1)*cos(2*pi*f(2)*temp_time_1st+A(2))+A(3)-temp_i_2nd;
    %x_2nd=lsqnonlin(fun_2nd,[-x0, p0, o0]);
    %x_2nd=lsqnonlin(fun_2nd,[-x0, p0, o0],[],[],[],[],[],[],@nonlcon_2nd);
    %background_lsq_half(nidx(i):nidx(i)+39)=x_1st(1)*cos(2*pi*f(2)*temp_time_1st+x_1st(2))+x_1st(3);
    %background_lsq_half(nidx(i)+40:nidx(i)+79)=x_2nd(1)*cos(2*pi*f(2)*temp_time_1st+x_2nd(2))+x_2nd(3);
    %faradic_lsq_half(nidx(i):nidx(i+1)-1)=current((nidx(i):nidx(i+1)-1))-background_lsq_half(nidx(i):nidx(i+1)-1);
end



%temp_i=current(nidx(3):nidx(3+1)-1);
%temp_time=time(nidx(3):nidx(3+1)-1)-time(nidx(3));

%fun = @(A)A(1)*cos(2*pi*f(2)*temp_time+A(2))+A(3)-temp_i;

%x0=abs(ft_i(2,1));
%p0=angle(ft_i(2,1));
%o0=mean(temp_i);

%x=lsqnonlin(fun,[x0, p0, o0]);

figure(100);
clf;
%plot(temp_time, temp_i,'k.-',temp_time, x(1)*cos(2*pi*f(2)*temp_time+x(2))+x(3),'b.-');
plot(time, current,'k.-',time, background_lsq,'b.-');

figure(101);
clf;
kk=cv_idx;
show_n=1;
plot(potential(nidx(kk):nidx(kk+show_n)-1), current(nidx(kk):nidx(kk+show_n)-1),'k.-');
%hold on;
%plot(potential(nidx(kk):nidx(kk+show_n)-1), background_lsq(nidx(kk):nidx(kk+show_n)-1),'b.-');
%hold on;
%plot(potential(nidx(kk):nidx(kk+show_n)-1), faradic_lsq(nidx(kk):nidx(kk+show_n)-1),'r.-');
axis([-0.2 0.6 -0.35 0.35]);
xticklabels({});
yticklabels({});
%xlabel('Potential') 
%ylabel('Current') 

figure(101);
clf;
kk=cv_idx;
show_n=1;
plot(abs(ft_i(:,kk)),'k.-');
axis([0 20 0 11]);
xticklabels({});
yticklabels({});
%xlabel('Potential') 
%ylabel('Current') 



figure(102);
clf;
%plot(temp_time, temp_i,'k.-',temp_time, x(1)*cos(2*pi*f(2)*temp_time+x(2))+x(3),'b.-');
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

%----------------------------------------
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
