time=fscv.time;
potential=fscv.potential;
current=fscv.current;


Fs = 1/(time(3)-time(2));
window_size=(idx(2)-idx(1));

nidx=iidx;
for i=2:length(idx)
    nidx(i)=window_size*(i-1)+iidx(1);
end


potential_fscv=potential(nidx(1):nidx(2)-1);
potential_fscv(42:end)=2*potential_fscv(41)-potential_fscv(42:end);

for i=1:length(nidx)-1
    current_fscv(:,i)=current(nidx(i):nidx(i+1)-1);
    faradic_fscv(:,i)=faradic_lsq(nidx(i):nidx(i+1)-1);
end

mag_df=1;

[ft_v,f,t] = stft(potential(iidx(1):iidx(end)),Fs, Window=rectwin(window_size), OverlapLength=0, FFTLength=window_size*mag_df, FrequencyRange="twosided");
[ft_i,f,t]=stft(current(iidx(1):iidx(end)),Fs, Window=rectwin(window_size), OverlapLength=0, FFTLength=window_size*1*mag_df, FrequencyRange="twosided");


temp_max=0;
temp_min=0;
temp1_max=0;
temp1_min=0;

for i=1:10
    temp_max=max(potential(iidx(i):idx(i+1)));
    temp_min=min(potential(iidx(i):iidx(i+1)));
    temp1_max=max(ifft(ft_v(:,i),'symmetric'));
    temp1_min=min(ifft(ft_v(:,i),'symmetric'));

end

global_scaler= (temp_max-temp_min)/(temp1_max-temp1_min);
global_offset=(temp_max+temp_min)/2-global_scaler*(temp1_max+temp1_min)/2;


%%
ft_z=ft_v./ft_i;
sdb = mag2db(abs(ft_z));
%sdb = abs(ft_z);
%sdb=angle(ft_z)*180/pi;
figure(10);
mesh(t,f,sdb);
%axis([0 100 0 500 -90 -10])
%cc = max(sdb(:))+[-35 10];
%cc = [0 20];
%cc = [-80 -70];
%cc = [250 275];
%cc = [0 0.5];
ax = gca;
ax.CLim = cc;
colormap jet;
view(2)
colorbar

figure(16);clf;plot(t,angle(ft_z(2,:))*180/pi);
figure(17);clf;plot(t,mag2db(abs(ft_z(2,:))));

f_idx=4;
figure(18);plot(t,angle(ft_i(f_idx,:))*180/pi);
%axis([0 100 -170 -120]);
figure(19);plot(t,mag2db(abs(ft_i(f_idx,:))));
%axis([0 100 -15 -1]);


% 특정시간 tt에 가장 가까운 CV를 찾아보자
tt=12.4;
[time_idx,cv_idx]=time2cv(tt,fscv,nidx);

tt1=26.404;
[time1_idx,cv1_idx]=time2cv(tt1,fscv,nidx);

tt2=27.74;
[time2_idx,cv2_idx]=time2cv(tt1,fscv,nidx);

%show_n=10;
%figure(22);clf;
%plot(potential(idx(cv_idx):idx(cv_idx+show_n)),faradic(idx(cv_idx):idx(cv_idx+show_n)),'r-');
%hold on;
%plot(fscv.potential(idx(cv_idx):idx(cv_idx+show_n)),fscv.current(idx(cv_idx):idx(cv_idx+show_n)),'k-');

show_nu=5;
%figure(23);clf;
%plot(potential(idx(cv1_idx):idx(cv1_idx+show_nu)),faradic(idx(cv1_idx):idx(cv1_idx+show_nu)),'r-');
%hold on;
%plot(potential(idx(cv1_idx):idx(cv1_idx+show_nu)),current(idx(cv1_idx):idx(cv1_idx+show_nu)),'k-');

figure(24);clf;
plot(fscv.time(idx(cv1_idx):idx(cv1_idx+show_nu)), fscv.current(idx(cv1_idx):idx(cv1_idx+show_nu)),'k-');
hold on;
plot(fscv.time(idx(cv1_idx):idx(cv1_idx+show_nu)), fscv.potential(idx(cv1_idx):idx(cv1_idx+show_nu))*0.27+0.05,'r-');


%sdb = abs(ft_z);
sdb1=angle(ft_z)*180/pi;
figure(25);
mesh(t,f,sdb1);

%cc = max(sdb(:))+[-35 10];
%cc = [2 20];

cc = [-90 -10];
%cc = [-90 -60];

%cc = [250 275];
%cc = [0 0.5];
ax = gca;
ax.CLim = cc;
view(2)
colorbar

%CV를 charging current에서 빼기 filter/substraction/guessing the background
%current

%global scaler and offset are values for ifft from fft.
temp_max=0;
temp_min=0;
temp1_max=0;
temp1_min=0;

for i=1:10
    temp_max=max(potential(idx(i):idx(i+1)));
    temp_min=min(potential(idx(i):idx(i+1)));
    temp1_max=max(ifft(ft_v(:,i),'symmetric'));
    temp1_min=min(ifft(ft_v(:,i),'symmetric'));

end

global_scaler= (temp_max-temp_min)/(temp1_max-temp1_min);
global_offset=(temp_max+temp_min)/2-global_scaler*(temp1_max+temp1_min)/2;


faradic=current;
background=current;

rmse_value=zeros(length(nidx)-1,1);
corr_value=zeros(length(nidx)-1,1);
mae_value=zeros(length(nidx)-1,1);

% without - in potential 초기버전
%local scaler : faradic containing one ferrocyanide
%36:40;;;36:44;;;;;;49:52;;;;;;;51:54;;;;56:57
%local scaler : nonfaradic no ferrocyanide 17:23,,,36:37

% with inversion
% 1:6;;;;;;;;;;

for i=1:length(nidx)-1
    temp_ft=ft_i(:,i);
    temp_ft(3:end)=0;
    temp_ift=ifft(temp_ft,'symmetric')*global_scaler+global_offset;
    temp_i=current(nidx(i):nidx(i+1)-1);
    local_scaler=mean(temp_i(36:37)./temp_ift(36:37));
    local_offset=mean(temp_i)-mean(temp_ift*local_scaler);
    background(nidx(i):nidx(i+1)-1)=local_scaler*temp_ift+local_offset;
    faradic(nidx(i):nidx(i+1)-1)=temp_i-local_scaler*temp_ift;
    rmse_value(i)=rmse(temp_i, local_scaler*temp_ift+local_offset);
    mae_value(i)=mae(temp_i, local_scaler*temp_ift+local_offset, ones(8,1));
    corr_value(i)=corr(temp_i,local_scaler*temp_ift+local_offset);
end

% half & half for cathodic anodic indepent local scaler

%for i=1:length(nidx)-1
%    temp_ft=ft_i(:,i);
%    temp_ft(3:end)=0;
%    temp_ift=ifft(temp_ft,'symmetric')*global_scaler+global_offset;
%    temp_i=current(nidx(i):nidx(i+1)-1);
%    local_scaler=mean(temp_i(36:37)./temp_ift(36:37));
%    local_offset=mean(temp_i)-mean(temp_ift*local_scaler);
%    background(nidx(i):nidx(i+1)-1)=local_scaler*temp_ift+local_offset;
%    faradic(nidx(i):nidx(i+1)-1)=temp_i-local_scaler*temp_ift;
%end

%%
figure(26);
clf;
kk=cv_idx;
show_n=1;
plot(potential(nidx(kk):nidx(kk+show_n)-1), current(nidx(kk):nidx(kk+show_n)-1),'k.-');
hold on;
plot(potential(nidx(kk):nidx(kk+show_n)-1), background(nidx(kk):nidx(kk+show_n)-1),'b.-');
hold on;
plot(potential(nidx(kk):nidx(kk+show_n)-1), faradic(nidx(kk):nidx(kk+show_n)-1),'r.-');
hold on;
plot(potential(nidx(kk):nidx(kk)+4), faradic(nidx(kk):nidx(kk)+4),'ro');

figure(27);
clf;
plot(time, faradic,'k.-');
axis([0 99.9999 -0.2 0.2]);

[ft_background,f,t] = stft(background(idx(1):idx(end)),Fs, Window=rectwin(window_size), OverlapLength=0, FFTLength=window_size*mag_df, FrequencyRange="onesided");
[ft_faradic,f,t]=stft(faradic(idx(1):idx(end)),Fs, Window=rectwin(window_size), OverlapLength=0, FFTLength=window_size*1*mag_df, FrequencyRange="onesided");

figure(28);
clf;
sdb1=angle(ft_background)*180/pi;
mesh(t,f,sdb1);
cc = [0 90];
ax = gca;
ax.CLim = cc;
view(2)
colorbar

figure(29);clf;plot(t,angle(ft_background(3,:))*180/pi);
figure(30);clf;plot(t,mag2db(abs(ft_background(3,:))));
figure(31);clf;plot(t,angle(ft_faradic(3,:))*180/pi);
figure(32);clf;plot(t,mag2db(abs(ft_faradic(3,:))));


%% -----------------simple substraction------------
faradic_sub=current;
temp_previous=current(nidx(1):nidx(2)-1);
for i=2:length(nidx)-1
    temp_i=current(nidx(i):nidx(i+1)-1);
    faradic_sub(nidx(i):nidx(i+1)-1)=temp_i-temp_previous;
    temp_previous=temp_i;
end

figure(30);
clf;
plot(time, faradic_sub,'k.-');
axis([0 49.9999 -0.1 0.1]);

figure(31);
clf;
kk=cv1_idx;
show_n=50;
plot(potential(idx(kk):idx(kk+show_n)-1), current(idx(kk):idx(kk+show_n)-1),'k.-');
hold on;
plot(potential(idx(kk):idx(kk+show_n)-1), faradic_sub(idx(kk):idx(kk+show_n)-1),'r.-');

figure(31);
mesh(t,potential_fscv,faradic_fscv,'LineStyle','none','FaceColor','interp');

cc = [-0.3 0.3];
ax = gca;
ax.CLim = cc;
view(2)
colormap("hsv");


function [tt_idx,cvt_idx] = time2cv(x,fscv,idx)
    [min_abs,t_idx]=min(abs(fscv.time-x));
    [min_abs,cv_idx]=min(abs(idx-t_idx));
    tt_idx=t_idx;
    cvt_idx=cv_idx;
end

