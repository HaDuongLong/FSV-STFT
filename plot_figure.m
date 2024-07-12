time2nidx=time;

for i=1:length(nidx)-1
    time2nidx(nidx(i):nidx(i+1))=linspace(i,i+1,nidx(2)-nidx(1)+1);
end

time2nidx(nidx(end):end)=linspace(time2nidx(nidx(end)),time2nidx(nidx(end))+1,length(time2nidx(nidx(end):end)));


z_resolution=0.001;
z_offset=1-round(min(faradic_lsq(nidx(1):nidx(end-1)))/z_resolution);
mesh_p=potential(nidx(1):nidx(2)-1);
mesh_t=time(nidx(1:end-1));
mesh_vol=ones(length(mesh_t), length(mesh_p)/2, int16((max(faradic_lsq(nidx(1):nidx(end-1)))-min(faradic_lsq(nidx(1):nidx(end-1))))/z_resolution), 'logical');

z_offset=1-round(min(current(nidx(1):nidx(end-1)))/z_resolution);
mesh_vol_raw=ones(length(mesh_t), length(mesh_p)/2, int16((max(current(nidx(1):nidx(end-1)))-min(current(nidx(1):nidx(end-1))))/z_resolution), 'logical');

for i=1:length(nidx)-1
    for k=1:40
        mesh_vol(i,k,round(faradic_lsq(nidx(i)+k-1)/z_resolution)+z_offset)=0;
        mesh_vol(i, 41-k,round(faradic_lsq(nidx(i)+k+39)/z_resolution)+z_offset)=0;
        mesh_vol_raw(i,k,round(current(nidx(i)+k-1)/z_resolution)+z_offset)=0;
        mesh_vol_raw(i, 41-k,round(current(nidx(i)+k+39)/z_resolution)+z_offset)=0;

    end
end


charge_anodic=time;
charge_cathodic=time;
charge_anodic_raw=time;
charge_back_anodic_raw=time;
charge_cathodic_raw=time;
charge_abs_raw=zeros(length(time),1);

for i=1:length(nidx)-1
    charge_anodic(nidx(i):nidx(i+1)-1)=sum(faradic_lsq(nidx(i):nidx(i)+38))*100*10^-4;
    charge_anodic_raw(nidx(i):nidx(i+1)-1)=sum(current(nidx(i):nidx(i)+38))*100*10^-4;
    charge_back_anodic_raw(nidx(i):nidx(i+1)-1)=sum(background_lsq(nidx(i):nidx(i)+39))*100*10^-4;
    charge_cathodic(nidx(i):nidx(i+1)-1)=sum(faradic_lsq(nidx(i)+39:nidx(i+1)-1))*100*10^-4;
    charge_cathodic_raw(nidx(i):nidx(i+1)-1)=sum(current(nidx(i)+39:nidx(i+1)-1))*100*10^-4;
    charge_abs_raw(nidx(i):nidx(i+1)-1)=(sum(background_lsq(nidx(i):nidx(i)+39))*100*10^-4-sum(current(nidx(i)+39:nidx(i+1)-1))*100*10^-4)/2;
end

charge_anodic(nidx(end):end)=0;
charge_anodic_raw(nidx(end):end)=0;
charge_back_anodic_raw(nidx(end):end)=0;
charge_cathodic(nidx(end):end)=0;
charge_cathodic_raw(nidx(end):end)=0;


faradic_lsq_anodic=zeros(length(nidx)-1,41);
faradic_lsq_cathodic=zeros(length(nidx)-1,41);
faradic_lsq_total=zeros(length(nidx)-1,82);
current_raw_anodic=zeros(length(nidx)-1,41);
current_raw_cathodic=zeros(length(nidx)-1,41);
current_raw_total=zeros(length(nidx)-1,82);
time_begin=zeros(length(nidx)-1,1);
potential_3d=potential(nidx(1):nidx(1)+40);
potential_3d_total=potential(nidx(1)-1:nidx(1)+80);

for i=1:length(nidx)-1
    faradic_lsq_anodic(i,:)=faradic_lsq(nidx(i):nidx(i)+40);
    faradic_lsq_cathodic(i,:)=faradic_lsq(nidx(i)+39:nidx(i+1)-1);
    faradic_lsq_total(i,:)=faradic_lsq(nidx(i)-1:nidx(i+1));
    current_raw_total(i,:)=current(nidx(i)-1:nidx(i+1));
    time_begin(i)=time(nidx(i));
end


%% charge -fitting전하 피팅
temp_x=time(nidx(3510:3700));
temp_y=charge_anodic_raw(nidx(3510:3700))*1000;
temp_x_plot=time(nidx(3510-20:3700+10));

temp_x1=time(nidx(3340:3470));
temp_y1=charge_anodic_raw(nidx(3340:3470))*1000;
temp_x1_plot=time(nidx(3340-25:3470+10));

temp_x2=time(nidx(3488:3494));
temp_y2=charge_anodic_raw(nidx(3488:3494))*1000;
temp_x2_plot=time(nidx(3488-3:3494+10))-temp_x2(1);

temp_x3=time(nidx(3317:3324));
temp_y3=charge_anodic_raw(nidx(3317:3324))*1000;
temp_x3_plot=time(nidx(3317-2:3324+10))-temp_x3(1);

figure(1000);clf;
plot(time,charge_anodic_raw*1000,'k.-','LineWidth',1,'MarkerSize',10);
hold on;
plot(temp_x_plot, 3.3284*10^8*exp(-0.6395.*temp_x_plot)+64.3090,'r-','LineWidth',1);
hold on;
plot(temp_x1_plot, 1.8001*10^8*exp(-0.6821.*temp_x1_plot)+63.8655,'r-','LineWidth',1);
hold on;
plot(temp_x2_plot+temp_x2(1), -0.3803*exp(-58.*temp_x2_plot)+69.6506,'g-','LineWidth',1);
hold on;
plot(temp_x3_plot+temp_x3(1), -0.3384*exp(-56.*temp_x3_plot)+66.1738,'g-','LineWidth',1);



figure(300);clf;
surf(time_begin(3350:3550,:),potential_3d,faradic_lsq_anodic(3350:3550,:)','EdgeColor','none','FaceColor','interp');
hold on;
surf(time_begin(3350:3550,:),potential_3d,faradic_lsq_cathodic(3350:3550,:)','EdgeColor','none','FaceColor','interp');
colormap jet;




figure(301);clf;
surf(time_begin(3400:3550,:),potential_3d_total,faradic_lsq_total(3400:3550,:)','EdgeColor','none','FaceColor','interp','FaceAlpha',0.8);
hold on;
plot3(time(nidx(3400)-1:nidx(3401)),potential_3d_total,faradic_lsq_total(3400,:)','k-','LineWidth',2);
hold on;
plot3(time(nidx(3483)-1:nidx(3484)),potential_3d_total,faradic_lsq_total(3483,:)','k-','LineWidth',2);
hold on;
plot3(time(nidx(3492)-1:nidx(3493)),potential_3d_total,faradic_lsq_total(3492,:)','k-','LineWidth',3);
colormap jet;
set(gca,'Ydir','reverse');

figure(302);clf;
mesh(time_begin(3400:3550,:),potential_3d_total,faradic_lsq_total(3400:3550,:)','EdgeColor','k','FaceColor','interp','LineStyle',':','FaceAlpha',0.8);
hold on;
plot3(time(nidx(3400)-1:nidx(3401)),potential_3d_total,faradic_lsq_total(3400,:)','k-','LineWidth',2);
hold on;
plot3(time(nidx(3483)-1:nidx(3484)),potential_3d_total,faradic_lsq_total(3483,:)','k-','LineWidth',2);
hold on;
plot3(time(nidx(3492)-1:nidx(3493)),potential_3d_total,faradic_lsq_total(3492,:)','k-','LineWidth',3);
hold on;
plot3(time(nidx(3280)-1:nidx(3281)),potential_3d_total,faradic_lsq_total(3280,:)','k-','LineWidth',2);
axis([27 28.4 -0.2 0.6 -0.12 0.12 ])
colormap jet;
set(gca,'Ydir','reverse');

figure(303);clf;
mesh(time_begin(3250:3550,:),potential_3d_total,faradic_lsq_total(3250:3550,:)'*100,'EdgeColor','k','FaceColor','interp','LineStyle',':','FaceAlpha',0.8);
hold on;
plot3(time(nidx(3250)-1:nidx(3251)),potential_3d_total,faradic_lsq_total(3250,:)'*100,'k-','LineWidth',3);
hold on;
plot3(time(nidx(3312)-1:nidx(3313)),potential_3d_total,faradic_lsq_total(3312,:)'*100,'b-','LineWidth',3);
hold on;
plot3(time(nidx(3323)-1:nidx(3324)),potential_3d_total,faradic_lsq_total(3323,:)'*100,'b-','LineWidth',3);
hold on;
plot3(time(nidx(3483)-1:nidx(3484)),potential_3d_total,faradic_lsq_total(3483,:)'*100,'r-','LineWidth',3);
hold on;
plot3(time(nidx(3492)-1:nidx(3493)),potential_3d_total,faradic_lsq_total(3492,:)'*100,'r-','LineWidth',3);
axis([25.9 28.4 -0.2 0.6 -12 12 ])
light("Style","local","Position",[20 100 0.5]);
colormap jet;
colorbar;
set(gca,'Ydir','reverse');

figure(304);clf;
mesh(time_begin(3250:3550,:),potential_3d_total,current_raw_total(3250:3550,:)'*100,'EdgeColor','k','FaceColor','interp','LineStyle',':','FaceAlpha',0.8);
hold on;
plot3(time(nidx(3250)-1:nidx(3251)),potential_3d_total,current_raw_total(3250,:)'*100,'k-','LineWidth',3);
hold on;
plot3(time(nidx(3312)-1:nidx(3313)),potential_3d_total,current_raw_total(3312,:)'*100,'b-','LineWidth',3);
hold on;
plot3(time(nidx(3323)-1:nidx(3324)),potential_3d_total,current_raw_total(3323,:)'*100,'b-','LineWidth',3);
hold on;
plot3(time(nidx(3483)-1:nidx(3484)),potential_3d_total,current_raw_total(3483,:)'*100,'r-','LineWidth',3);
hold on;
plot3(time(nidx(3492)-1:nidx(3493)),potential_3d_total,current_raw_total(3492,:)'*100,'r-','LineWidth',3);
axis([25.9 28.4 -0.2 0.6 -40 40 ])
light("Style","local","Position",[20 100 0.5]);
colormap jet;
colorbar;
set(gca,'Ydir','reverse');


figure(305);clf;
plot(potential_3d_total,faradic_lsq_total(3311,:)'*100,'k--','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3314,:)'*100,'r--','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3315,:)'*100,'b--','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3317,:)'*100,'k-','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3319,:)'*100,'r-','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3323,:)'*100,'b-','LineWidth',1);
hold on;
axis([-0.2 0.6 -8 8 ])


figure(306);clf;
plot(potential_3d_total,faradic_lsq_total(3482,:)'*100,'k--','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3483,:)'*100,'r--','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3484,:)'*100,'b--','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3485,:)'*100,'k-','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3486,:)'*100,'r-','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(3492,:)'*100,'b-','LineWidth',1);
hold on;
axis([-0.2 0.6 -12 12 ]);

figure(307);clf;
surf(time_begin(3250:3550,:),potential_3d_total,faradic_lsq_total(3250:3550,:)'*100,'EdgeColor','none','FaceColor','interp');
hold on;
plot3(time(nidx(3250)-1:nidx(3251)),potential_3d_total,faradic_lsq_total(3250,:)'*100,'k-','LineWidth',3);
hold on;
plot3(time(nidx(3312)-1:nidx(3313)),potential_3d_total,faradic_lsq_total(3312,:)'*100,'b-','LineWidth',3);
hold on;
plot3(time(nidx(3323)-1:nidx(3324)),potential_3d_total,faradic_lsq_total(3323,:)'*100,'b-','LineWidth',3);
hold on;
plot3(time(nidx(3483)-1:nidx(3484)),potential_3d_total,faradic_lsq_total(3483,:)'*100,'r-','LineWidth',3);
hold on;
plot3(time(nidx(3492)-1:nidx(3493)),potential_3d_total,faradic_lsq_total(3492,:)'*100,'r-','LineWidth',3);
axis([25.9 28.4 -0.2 0.6 -12 12 ])
colormap jet;
colorbar;
view(90,0);
set(gca,'Ydir','reverse');


figure(402);clf;
mesh(time_begin(1250:1450,:),potential_3d_total,current_raw_total(1250:1450,:)'*100,'EdgeColor','k','FaceColor','interp','LineStyle',':','FaceAlpha',0.8);
hold on;
plot3(time(nidx(1250)-1:nidx(1251)),potential_3d_total,current_raw_total(1250,:)'*100,'k-','LineWidth',2);
hold on;
plot3(time(nidx(1330)-1:nidx(1331)),potential_3d_total,current_raw_total(1330,:)'*100,'k-','LineWidth',2);
axis([9.99 11 -0.2 0.6 -15 15])
colormap jet;
colorbar;
set(gca,'Ydir','reverse');

figure(403);clf;
surf(time_begin(1250:1450,:),potential_3d_total,faradic_lsq_total(1250:1450,:)'*100,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.7);
hold on;
plot3(time(nidx(1250)-1:nidx(1251)),potential_3d_total,faradic_lsq_total(1250,:)'*100,'k-','LineWidth',2);
hold on;
plot3(time(nidx(1330)-1:nidx(1331)),potential_3d_total,faradic_lsq_total(1330,:)'*100,'k-','LineWidth',2);
axis([9.99 11 -0.2 0.6 -10 2])
colormap jet;
colorbar;
set(gca,'Ydir','reverse');

figure(404);clf;
plot(potential_3d_total,faradic_lsq_total(1330,:)'*100,'r-','LineWidth',1.5);
hold on;
plot(potential_3d_total,current_raw_total(1330,:)'*100,'b-','LineWidth',1);
hold on;
plot(potential_3d_total,current_raw_total(1300,:)'*100,'k--','LineWidth',1);
hold on;
plot(potential_3d_total,faradic_lsq_total(1300,:)'*100,'k.-','LineWidth',1);

figure(405);clf;
mesh(time_begin(550:700,:),potential_3d_total,current_raw_total(550:700,:)'*100,'EdgeColor','k','FaceColor','interp','LineStyle',':','FaceAlpha',0.8);
hold on;
plot3(time(nidx(550)-1:nidx(551)),potential_3d_total,current_raw_total(550,:)'*100,'k-','LineWidth',2);
hold on;
plot3(time(nidx(616)-1:nidx(617)),potential_3d_total,current_raw_total(616,:)'*100,'k-','LineWidth',2);
colormap jet;
colorbar;
set(gca,'Ydir','reverse');

figure(406);clf;
surf(time_begin(550:700,:),potential_3d_total,faradic_lsq_total(550:700,:)'*100,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.7);
hold on;
plot3(time(nidx(550)-1:nidx(551)),potential_3d_total,faradic_lsq_total(550,:)'*100,'k-','LineWidth',2);
hold on;
plot3(time(nidx(616)-1:nidx(617)),potential_3d_total,faradic_lsq_total(616,:)'*100,'k-','LineWidth',2);
colormap jet;
colorbar;
set(gca,'Ydir','reverse');


figure(200);
clf;
plot(time2nidx,current,'k.-');


figure(201);
clf;
iii=2:100:500;
for i=1:length(iii)
    plot3(time(nidx(iii(i)):nidx(iii(i)+1)), potential(nidx(iii(i)):nidx(iii(i)+1)), faradic_lsq(nidx(iii(i)):nidx(iii(i)+1))*100,'k.-');
    axis([0 80 -0.2 0.7 -10 10]);
    hold on;
end

iii=550:300:1570;
for i=1:length(iii)
    plot3(time(nidx(iii(i)):nidx(iii(i)+1)), potential(nidx(iii(i)):nidx(iii(i)+1)), faradic_lsq(nidx(iii(i)):nidx(iii(i)+1))*100,'y.-');
    axis([0 80 -0.2 0.7 -10 10]);
    hold on;
end

iii=1570:2:1575;
for i=1:length(iii)
    plot3(time(nidx(iii(i)):nidx(iii(i)+1)), potential(nidx(iii(i)):nidx(iii(i)+1)), faradic_lsq(nidx(iii(i)):nidx(iii(i)+1))*100,'b.-');
    axis([0 80 -0.2 0.7 -10 10]);
    hold on;
end


iii=3310:30:4000; 
for i=1:length(iii)
    plot3(time(nidx(iii(i)):nidx(iii(i)+1)), potential(nidx(iii(i)):nidx(iii(i)+1)), faradic_lsq(nidx(iii(i)):nidx(iii(i)+1))*100,'r.-');
    axis([0 80 -0.2 0.7 -15 15]);
    hold on;
end



figure(202);
clf;
iii=2:100:5000;
for i=1:length(iii)
    plot3(time(nidx(iii(i)):nidx(iii(i)+1)), potential(nidx(iii(i)):nidx(iii(i)+1)), faradic_lsq(nidx(iii(i)):nidx(iii(i)+1))*100);
    axis([0 80 -0.2 0.7 -10 10]);
    hold on;
end

figure(203);
clf;
iii=1200:1:1400;
for i=1:length(iii)
    plot3(time(nidx(iii(i)):nidx(iii(i)+1)), potential(nidx(iii(i)):nidx(iii(i)+1)), current(nidx(iii(i)):nidx(iii(i)+1))*100);
end
