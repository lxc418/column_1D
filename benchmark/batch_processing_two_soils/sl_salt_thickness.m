% clear
% load plot.mat
% fclose('all');
run('/storage/macondo/s4524462/SutraLab/mfiles/slsetpath.m')
c=ConstantObj();

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



    solidmass_matrix_kg = reshape(nod(time_step).terms{sm_idx},[inp.nn1,inp.nn2]);
    solidmass_surface_kg(1:inp.nn2) = solidmass_matrix_kg(inp.nn1,:);
    solidmass_thickness_mm(1:inp.nn2) = solidmass_surface_kg./c.density_solid_nacl_kgPm3./area1_m2*c.m2mm;

% %calculate average solid salt thickness	
left_solid_salt  = mean( solidmass_thickness_mm(1:(inp.nn2-1)/2) );
right_solid_salt = mean( solidmass_thickness_mm((inp.nn2-1)/2+1:end) );

average_solid_salt = [left_solid_salt right_solid_salt];
writematrix(average_solid_salt,'../M.xlsx','Sheet',1,'Range','aFINDMEROW:bFINDMEROW')

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

% figure
% [csal,hcsal] = contourf(x_matrix,y_matrix,c_matrix);
    % hold on
    % set(hcsal,'EdgeColor','none');
    % color = jet;
    % colormap(gca,color);
    % cbsal = colorbar;
	% caxis([0 0.246])
    % cbsal.Label.String = 'Concentration (-)';
% % scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
% savefig('concentration_contour.fig')						 

% figure
% [csal,hcsal] = contourf(x_matrix,y_matrix,s_matrix);
    % hold on
    % set(hcsal,'EdgeColor','none');
    % color = jet;
    % color = flipud(color);
    % colormap(gca,color);
    % cbsal = colorbar;
	% caxis([0 1])
    % cbsal.Label.String = 'Saturation (-)';
% % scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
% savefig('saturation_contour.fig')	
