clear
load plot.mat
fclose('all');
c=ConstantObj();

time_step = length(bcof);
time_day  = [bcof.tout]/3600/24;%second to day
time_nod_day = arrayfun(@(y) y.tout,nod) * c.dayPsec;
water_table  = inp.pbc/9800;

x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);%inp.nn2 is number of nodes in x direction 
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);
x_ele_matrix = reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
y_ele_matrix = reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);

%% evaporation data
evapo_kgs = zeros(inp.nn2,time_step);
for i=1:inp.nn2
    
if i<inp.nn2   
    area1(1:i)    = (x_matrix(1,i+1)-x_matrix(1,i))*inp.z(1); %evaporation area 
else
    area1(1:i)    = (x_matrix(1,i)-x_matrix(1,i-1))*inp.z(1); %the right end node
end 

    evapo_kgs(i,:)  = -arrayfun(@(y) y.qin(i),bcof);
    evapo_mmday     = evapo_kgs/area1(i)*c.secPday; %evaporation rate of every surface node
    
end

evapo_mmday(1,:)    =   evapo_mmday(1,:)*2;
evapo_mmday(end,:)  =   evapo_mmday(end,:)*2;%avoid boundary effect by double the two nodes
total_evapo_mmday(1:time_step)   =  sum (evapo_mmday(:,1:time_step))./(length(area1)); %total evaporation rate

cumulative_evapo_mm =   zeros(1,time_step); %cumulative evaporation amount
for i=2:time_step
    cumulative_evapo_mm(1) =  total_evapo_mmday(1)*inp.scalt*inp.nbcfpr*c.dayPsec;
    cumulative_evapo_mm(i) =  total_evapo_mmday(i)*inp.scalt*inp.nbcfpr*c.dayPsec + cumulative_evapo_mm(i-1);
end

   

%% plot control
fig_pos.left   = 0.05;
fig_pos.bottom = 0.7;
fig_pos.length = 0.35;
fig_pos.height = 0.26;

%nt=10;
a.fs = 15;
a.lw = 2; %line width
a.cz = 8; %the size of the marker
fs   = 2; % sampling frequency
% Creating the movie : quality = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt = 100;
%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig = figure;
set (gcf,'Position',[0,0,1920,1080]); %resolution 1080p
% set(gcf,'Units','normalized', 'OuterPosition',[0 0 1 1]);  % maximize the plotting figure

mov           =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
mov.FrameRate = 5;
mov.Quality   = qt;
open(mov);

for nt=1:round(time_step/40):time_step
    %--- plot concentration---- 
    subplot('Position',[0.05 0.58 0.14 0.38]),
%     plot(clab(2,:,1),clab(1,:,1),'rd','MarkerSize',cz,'MarkerFaceColor','r'); hold on
%     plot(clab(2,:,2),clab(1,:,2),'go','MarkerSize',cz,'MarkerFaceColor','g');hold on
%     plot(clab(2,:,3),clab(1,:,3),'ks','MarkerSize',cz,'MarkerFaceColor','k');hold on;
%     hleg1 = legend('74% Sat.','50% Sat.','32% Sat.','Location','SouthEast');
%     set(hleg1, 'Box', 'on','FontSize',a.fz,'LineWidth',a.lw)
	c_matrix  = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]);
    plot(c_matrix(2,:),y_matrix(2,:),'LineWidth',a.lw);
    xlabel('Concentration (kg/kg)','FontSize',a.fz,'FontWeight','bold')
    ax1 = gca;
    set(ax1,'FontSize',fl,'FontWeight','bold','LineWidth',a.lw)
    ylabel('Depth (m)','FontSize',a.fz,'FontWeight','bold')
    axis([0 0.3 y_matrix(2,1) y_matrix(2,end)])
    %% ---plot Saturation-----
    subplot('Position',[0.20 0.58 0.14 0.38]),
%     plot(slab(2,:,1),slab(1,:,1),'rd','MarkerSize',cz,'MarkerFaceColor','r');hold on
%     plot(slab(2,:,2),slab(1,:,2),'go','MarkerSize',cz,'MarkerFaceColor','g');hold on
%     plot(slab(2,:,3),slab(1,:,3),'kS','MarkerSize',cz,'MarkerFaceColor','k');hold on
	s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);
    plot(s_matrix(2,:),y_matrix(2,:),'LineWidth',a.lw);
    xlabel('Saturation (m^{3} m^{-3})','FontSize',a.fz,'FontWeight','bold')
    ax1 = gca;
    set(ax1,'YTickLabel','','FontSize',fl,'FontWeight','bold','LineWidth',a.lw)
%     hleg1 = legend('74% Sat.','50% Sat.','32% Sat.','Location','SouthWest');
%     set(hleg1, 'Box', 'on','FontSize',a.fz,'LineWidth',a.lw)
    axis([0 1 y_matrix(2,1) y_matrix(2,end)])    
   %% ---plot temperature-----
    subplot('Position',[0.35 0.58 0.14 0.38]),
%     plot(tlab(2,:,1),tlab(1,:,1),'rd','MarkerSize',cz,'MarkerFaceColor','r');hold on
%     plot(tlab(2,:,2),tlab(1,:,2),'go','MarkerSize',cz,'MarkerFaceColor','g');hold on
%     plot(tlab(2,:,3),tlab(1,:,3),'ks','MarkerSize',cz,'MarkerFaceColor','k');hold on
    temp_matrix = reshape(nod(nt).terms{temp_idx},[inp.nn1,inp.nn2])-273.15;    
    plot(temp_matrix(2,:),y_matrix(2,:),'LineWidth',a.lw);
    xlabel('Temperature (\circC)','FontSize',a.fz,'FontWeight','bold')
    axis([20 50 y_matrix(2,1) y_matrix(2,end)])    
    ax1 = gca;
    set(ax1,'YTickLabel','','FontSize',fl,'FontWeight','bold','LineWidth',a.lw)    
%     hleg1 = legend('74% Sat.','50% Sat.','32% Sat.','Location','SouthEast');
%     set(hleg1, 'Box', 'on','FontSize',a.fz,'LineWidth',a.lw)
    
    %% ----Total Head output-------
    subplot('Position',[0.50 0.58 0.14 0.38])
	p_matrix  = reshape(nod(nt).terms{p_idx},[inp.nn1,inp.nn2]);
	total_head_matrix=p_matrix(2,:)/c.g./(1000+702.24*(c_matrix(2,:)-0.0001))+y_matrix(2,:);
    plot(total_head_matrix,y_matrix(2,:),'LineWidth',a.lw);
    xlabel('Total Head (m)','FontSize',a.fz,'FontWeight','bold')
    axis([-5000 0 y_matrix(2,1) y_matrix(2,end)])
    ax1 = gca;
    set(ax1,'YTickLabel','','FontSize',fl,'FontWeight','bold','LineWidth',a.lw)
    %set(ax1,'YAxisLocation','right','XAxisLocation','top')
    
    %% ----Vertical flux output-------
    subplot('Position',[0.65 0.58 0.14 0.38])
	vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
    plot(vy_matrix(2,:)*c.secPday*c.m2mm,y_ele_matrix(2,:),'r',qv.qvy(2,:)*c.secPday*c.m2mm,y_ele_matrix(2,:)),'b','LineWidth',a.lw) 
    xlabel('Vertical Flux (mm day^{-1})','FontSize',a.fz,'FontWeight','bold')
    hleg1 = legend('liquid water','water vapor','Location','SouthEast');
    set(hleg1, 'Box', 'on','FontSize',a.fz)
    axis([-1 5 y_ele_matrix(2,1) y_ele_matrix(2,end)])     
    ax1 = gca;
    set(ax1,'YTickLabel','','FontSize',fl,'FontWeight','bold','LineWidth',a.lw)
    %% ----Relative K output-------
    subplot('Position',[0.80 0.58 0.14 0.38])
	kr_matrix = reshape(ele(nt).terms{kr_idx},[inp.nn1-1,inp.nn2-1]);
    semilogx(kr_matrix(2,:), y_ele_matrix(2,:)),'LineWidth',a.lw)
    xlabel('Relative K','FontSize',a.fz,'FontWeight','bold')
    ylabel('Depth (m)','FontSize',a.fz,'FontWeight','bold')
    axis([10^-20 1 y_ele_matrix(2,1) y_ele_matrix(2,end)])
    ax1 = gca;
    set(ax1,'YAxisLocation','right','YColor','k','FontSize',fl,'FontWeight','bold','LineWidth',a.lw)
    
    
    
	F = getframe(gcf); % save the current figure
    writeVideo(mov,F);% add it as the next frame of the movie

end
saveas(a.fig,'a','fig')
close(mov);
close(a.fig);
%figure
% [csal,hcsal] = contourf(x_matrix,y_matrix,p_mtx);
    % hold on
    % set(hcsal,'EdgeColor','none');
    % color = jet;
    % colormap(gca,color);
    % cbsal = colorbar;
    % cbsal.Label.String = 'Pressure (-)';
% scatter (a1(1,(1:f3(1)),1),a1(2,(1:f3(1)),1),2,'filled','w')
% savefig('pressure_contour.fig')						 

figure
[csal,hcsal] = contourf(x_matrix,y_matrix,c_matrix);
    hold on
    set(hcsal,'EdgeColor','none');
    color = jet;
    colormap(gca,color);
    cbsal = colorbar;
	caxis([0 0.246])
    cbsal.Label.String = 'Concentration (-)';
scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
savefig('concentration_contour.fig')						 

figure
[csal,hcsal] = contourf(x_matrix,y_matrix,s_matrix);
    hold on
    set(hcsal,'EdgeColor','none');
    color = jet;
    color = flipud(color);
    colormap(gca,color);
    cbsal = colorbar;
	caxis([0 1])
    cbsal.Label.String = 'Saturation (-)';
scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
savefig('saturation_contour.fig')	
