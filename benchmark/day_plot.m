%% this file is used for plot transient results at a given time(day) 
run('/storage/macondo/s4524462/SutraLab/mfiles/slsetpath.m')
c=ConstantObj();
day_output = [7];%set days need to be plotted	

time_step = length(et);
time_day  = [bcof.tout]/3600/24;%second to day
time_nod_day = arrayfun(@(y) y.tout,nod) * c.dayPsec;
water_table  = inp.pbc/(c.rhow_pure_water+700*0.035);	

x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);%inp.nn2 is number of nodes in y direction 
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);
x_ele_matrix = reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
y_ele_matrix = reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);

%locate the center of left & right for different soil types
[numRows,numCols] = size (x_matrix);
left_centre       = round((numCols+1)/4);
right_centre      = round((numCols+1)/4+(numCols-1)/2);
%% evaporation data (from bcof without the vapor contribution)
% evapo_kgs = zeros(time_step,inp.nn2);
			   
% for i=1:inp.nn2
																	   
% if i<inp.nn2   
    % area1_m2(1:i)    = (x_matrix(1,i+1)-x_matrix(1,i))*inp.z(1); %evaporation area 
% else
    % area1_m2(1:i)    = (x_matrix(1,i)-x_matrix(1,i-1))*inp.z(1); %the right end node
% end 

    % evapo_kgs(i,:)  = -arrayfun(@(y) y.qin(i),bcof);
    % evapo_mmday     = evapo_kgs/area1_m2(i)*86400; %evaporation rate of every surface node
    										   
% end

%% solute inflow from bottom (from bcof without the vapor contribution)

for i= 1:inp.nn2

    solute_kgs(i,:)  = -arrayfun(@(y) y.qpu(i),bcop);
    
end
solute_gday= solute_kgs'.*c.kg2g*c.secPday;

area1_m2    = (x_matrix(1,2)-x_matrix(1,1))*inp.z(1);
for i=2:time_step

	evapo_mmday(1,:)  		 =  reshape(et(1).terms{et_idx},[1,inp.nn2])*c.ms2mmday;
	total_evapo_mmday(1)     =  sum (evapo_mmday(1,:))./inp.nn2; %the evp rate from the whole surface
	evapo_mmday(i,:)  		 =  reshape(et(i).terms{et_idx},[1,inp.nn2])*c.ms2mmday;
	total_evapo_mmday(i)     =  sum (evapo_mmday(i,:))./inp.nn2;
	
    cumulative_evapo_mm(1)   =  total_evapo_mmday(1)*inp.scalt*inp.nbcfpr*c.dayPsec;
    cumulative_evapo_mm(i)   =  total_evapo_mmday(i)*inp.scalt*inp.nbcfpr*c.dayPsec + cumulative_evapo_mm(i-1);
   
end

   
%% plot control
fig_pos.left   = 0.1;
fig_pos.bottom = 0.755;
fig_pos.length = 0.07;
fig_pos.height = 0.2;

%nt=10;
a.fs = 15;
a.lw = 2; %line width
a.cz = 8; %the size of the marker
fs   = 5; % sampling frequency
% Creating the movie : quality = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt = 100;
%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig = figure;
set (gcf,'Position',[0,0,1920,1080]); %resolution 1080p
% set(gcf,'Units','normalized', 'OuterPosition',[0 0 1 1]);  % maximize the plotting figure

%calculate time_step for output
timestep_output = round(day_output*c.secPday/inp.nprint/inp.scalt)+1;%jump the first timestep which is 1 second
timestep_output (timestep_output>time_step)=time_step;
% timestep_output = time_step; %just plot the end of simulation

for nt = timestep_output
    s_matrix  = reshape(nod(nt+1).terms{s_idx},[inp.nn1,inp.nn2]);	 %write in matrix form.
    s_surface_matrix = s_matrix(inp.nn1,:);	   
	c_matrix = reshape(nod(nt+1).terms{c_idx},[inp.nn1,inp.nn2])*1000;%change unit to ppt
	c_surface_matrix = c_matrix(inp.nn1,:);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);

%% -------------  sub 1 Salinity & solid salt  ---------------------
    a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom+fig_pos.height/2+0.003,...
          fig_pos.length,fig_pos.height/2+0.03]);
yyaxis left	  
    a.plot1=plot(x_matrix(1,:), c_surface_matrix,...
             '-','linewidth',a.lw,'color',[0.4940 0.1840 0.5560]);hold off
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
	axis([0 x_matrix(1,end) 0 300])
    yticks([100,200,300])
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
	set(gca,'YColor',[0.4940 0.1840 0.5560]);
	ylabel({'Salinity (ppt)'},'FontSize',a.fs);
yyaxis right
    solidmass_matrix_kg = reshape(nod(nt+1).terms{sm_idx},[inp.nn1,inp.nn2]);
    solidmass_surface_kg(1:inp.nn2) = solidmass_matrix_kg(inp.nn1,:);
    solidmass_thickness_mm(1:inp.nn2) = solidmass_surface_kg./c.density_solid_nacl_kgPm3./area1_m2*c.m2mm;
    a.plot2 = plot(x_matrix(1,:), solidmass_thickness_mm(1:inp.nn2),...
             'r-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
    % ylabel({'Solid salt';'(mm)'},'FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 0.6])
    yticks([0,0.2,0.4,0.6])
	set(gca,'YColor','r');
	grid on
	xticks([0.05,0.1,0.15]);
	xticklabels([]);
%% -------------  sub 2 salinity contour  ---------------------
a.sub2=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-fig_pos.height/2,...
          fig_pos.length+0.0305,fig_pos.height]);
    
    % a.plot4=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold off	
	a.plot2=contourf(x_matrix,y_matrix,log10(c_matrix+1-inp.ubc(1)*1000),'LevelStep',0.005,'LineStyle','None');hold on
    color = jet;
    colormap(gca,color);
	caxis([0.09 2.35])
	set(gca,'ColorScale','log')%log scale of colormap
	cbsal = colorbar;
    cbsal.Label.String = 'Salinity (ppt)';
	cbsal.TicksMode    = 'manual';
	cbsal.Ticks        = [0.09,0.3,0.78,1.2,1.82,2.35];
	cbsal.TickLabels   = [35,36,40,50,100,260];	
	set(colorbar,'visible','off')
%	scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    % color = jet;
    % colormap(gca,color);
	% caxis([0 0.264])
    % cbsal = colorbar;
    % cbsal.Label.String = 'Concentration (-)';
    % get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    get(gca,'xtick');
	xticks([0.05,0.1,0.15]);
	xticklabels([]);
    % xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
	yticks([0.1,0.2,0.3,0.4])
	axis([0 x_matrix(1,end) 0 y_matrix(end,1)])       
	
 %% -------------  sub 3 sat & evt  ---------------------
    a.sub3=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-fig_pos.height-0.035,...
          fig_pos.length,fig_pos.height/2+0.03]);
yyaxis left
    a.plot3=plot(x_matrix(1,:), s_surface_matrix,...
             '-','linewidth',a.lw,'color',[0.4660 0.6740 0.1880]);hold off		  			 
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
	axis([0 x_matrix(1,end) 0 1])
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
	set(gca,'YColor',[0.4660 0.6740 0.1880]);
    yticks([0.5,1])
	ylabel({'Saturation(-)'},'FontSize',a.fs);		  
yyaxis right		   
    a.plot3=plot(x_matrix(1,:), evapo_mmday(nt,:),...
             '-','linewidth',a.lw,'color',[0 0.4470 0.7410]);hold off
	grid on
	% grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
    % xlabel('x','FontSize',a.fs);
	set(gca,'YColor',[0 0.4470 0.7410]);
	xticks([0.05,0.1,0.15]);
	xticklabels([]);
	yticks([0,8,16])
    % ylabel({'Evt';'(mm/day)'},'FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 16])  
%% -------- sub 4 saturation contour  ---------
a.sub4=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-0.037-fig_pos.height*2,...
          fig_pos.length+0.0305,fig_pos.height]);
    a.plot4=contourf(x_matrix,y_matrix,s_matrix,'EdgeColor','none');hold off

%    scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = hot;
    color = flipud(color);
    colormap(gca,color);
	caxis([0 1])
    cbsal = colorbar;
    cbsal.Label.String = 'Saturation (-)';
	set(colorbar,'visible','off')
    get(gca,'xtick');
	xticks([]);
    set(gca,'fontsize',a.fs);
	yticks([0.1,0.2,0.3,0.4])
    ylabel('Elevation (m)','FontSize',a.fs);
	axis([0 x_matrix(1,end) 0 inf])        

%% -------- sub 4 vapor & liquid flow ---------
a.sub4=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-0.04-fig_pos.height*3,...
          fig_pos.length,fig_pos.height]);		  
	a.plot4= quiver(x_ele_matrix(1:4:end,1:4:end),y_ele_matrix(1:4:end,1:4:end),vx_matrix(1:4:end,1:4:end),vy_matrix(1:4:end,1:4:end));%plot quiver with larger space
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
	hold on
    qvx_mtx = qv(nt).qvx; % vector of vapor is acquired from the concentration difference bewteen two nodes,
    qvy_mtx = qv(nt).qvy; % which is used to calculate the vector in element center.
	qvx_plot_mtx = zeros(inp.nn1-1,inp.nn2-1);  
    qvy_plot_mtx = zeros(inp.nn1-1,inp.nn2-1);
    for i = 1:inp.nn2-1 %along x
        for j = 1:inp.nn1-1 %along y
            qvx_plot_mtx(j,i) = (qvx_mtx(j,i)+qvx_mtx(j+1,i))/2;
            qvy_plot_mtx(j,i) = (qvy_mtx(j,i)+qvy_mtx(j,i+1))/2;
        end
    end
    % a.plot8=quiver(x_ele_matrix,y_ele_matrix,qvx_plot_mtx,qvy_plot_mtx,'k-');hold off
	a.plot4= quiver(x_ele_matrix(1:4:end,1:4:end),y_ele_matrix(1:4:end,1:4:end),qvx_plot_mtx(1:4:end,1:4:end),qvy_plot_mtx(1:4:end,1:4:end),'r');hold off%plot quiver with larger space
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');	
	xticks([0,0.05,0.1,0.15]);
    set(gca,'fontsize',a.fs);
	xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 y_matrix(inp.nn1,1)])	 
	
% %% -------- sub 5 vapor & liquid flow ---------
% a.sub4=subplot('position'...
         % ,[fig_pos.left,fig_pos.bottom-fig_pos.height*4,...
          % fig_pos.length,fig_pos.height]);		
		  
    % kr_matrix = reshape(ele(nt).terms{kr_idx},[inp.nn1-1,inp.nn2-1]);    
    % a.plot7=semilogx(kr_matrix(:,left_centre),y_ele_matrix(:,1),'-','color',[0 0.4470 0.7410],'linewidth',a.lw);hold on
    % a.plot7=semilogx(kr_matrix(:,right_centre),y_ele_matrix(:,1),'-','color',[0.8500 0.3250 0.0980],'linewidth',a.lw);hold off
    % get(gca,'xtick');
    % set(gca,'fontsize',a.fs);
    % xlabel('Relative K (-)','FontSize',a.fs);
    % ylabel('Elevation (m)','FontSize',a.fs);
    % axis([-inf 1.09 0 y_matrix(inp.nn1,1)])			  
		  
		  		 				  	  
end																

	figure_name=sprintf('day_%.2f.fig',nod(nt+1).tout*c.dayPsec);
	saveas(a.fig,figure_name)




