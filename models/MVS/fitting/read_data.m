% Written by Yan Peng, UBC, 2016/10/19
% Modified by Justin Reiher, UBC, 2017/05/30

function [idvd_n, idvd_p, idvg_n, idvg_p] = read_data(idvd_n_file, idvd_p_file, idvg_n_file, idvg_p_file, scale)
    
	%% Constants for Ids vs Vds plots
	% Start,stop and Vds_step size for Vds in Ids vs Vds curves
     Vds_start = scale(1);
     Vds_step = scale(2);
     Vds_ending = scale(3);
     Vds_length = length(Vds_start:Vds_step:Vds_ending);

	%Start, stop and Vgs_step size for number of Vgs steps in Ids vs Vds curves
     Vgs_start = scale(4);
     Vgs_step = scale(5);
     Vgs_ending = scale(6);
     num_Vgs_steps = length(Vgs_start:Vgs_step:Vgs_ending);

 
	% read in Ids vs Vds data for different Vgs for the nmos, need to add in Vgs values
	% form of the data is Vds, Vgs, Ids
    idvd_n_raw=dlmread(idvd_n_file);
    idvd_n = zeros(Vds_length*num_Vgs_steps,3);
    idvd_n(:,1) = idvd_n_raw(:,1);
    idvd_n(:,3) = idvd_n_raw(:,2);
    for j = 0:(num_Vgs_steps-1)
		for i = 1:Vds_length
			idvd_n(i+j*Vds_length,2) = Vgs_start + j*Vgs_step;
		end
    end


    %% read in Ids vs Vds data for different Vgs for the pmos
	%form of the data is Vds, Vgs, Ids
    idvd_p_raw=dlmread(idvd_p_file);
    idvd_p = zeros(Vds_length*num_Vgs_steps,3);

    %need to make Vds negative, and add in Vgs (also negative)
    %Data in form of Vds, Vgs, Ids
    idvd_p(:,1) = Vgs_ending-idvd_p_raw(:,1);
    idvd_p(:,2) = Vgs_ending-idvd_p_raw(:,4);
    idvd_p(:,3) = idvd_p_raw(:,6);

%     %%for j = 0:(num_Vgs_steps-1)
% 		for i = 1:Vds_length
% 			idvd_p(i+j*Vds_length,2) = -Vgs_start - j*Vgs_step;
% 		end
%     end	


   %% Constants for Ids vs Vgs plots
   %   
   Vds_step_start = scale(7);
   Vds_step_end = scale(8);

   Vgs_start = scale(9);
   Vgs_step_size = scale(10);
   Vgs_end = scale(11);
   Vgs_length = length(Vgs_start:Vgs_step_size:Vgs_end);

    %% Data for Ids vs Vgs for nmos
    % data in the form of Vgs, Ids @ Vds_step_star, Ids  @Vds_step_end
    idvg_n_raw=dlmread(idvg_n_file);
    idvg_n = zeros(Vgs_length,3);
    Vgs_length = length(Vgs_start:Vgs_step_size:Vgs_end);

    idvg_n(:,1) = Vgs_start:Vgs_step_size:Vgs_end;

    for j = 0:1
		for i = 1:Vgs_length
			idvg_n(i,j+2) = idvg_n_raw(i+j*Vgs_length,2);
		end
    end

    
    %% Data for Ids vs Vgs for pmos
    % data in the form of Vgs, Ids @ Vds_step_star, Ids  @Vds_step_end
    idvg_p_raw=dlmread(idvg_p_file);
    

    idvg_p = zeros(Vgs_length,3);
    Vgs_length = length(Vgs_start:Vgs_step_size:Vgs_end);

    idvg_p(:,1) = Vgs_start:Vgs_step_size:Vgs_end;
    idvg_p(:,1) = Vgs_end-idvg_p(:,1);

    for j = 0:1
		for i = 1:Vgs_length
			idvg_p(i,3-j) = idvg_p_raw(i+j*Vgs_length,6);
		end
    end

idvg_p(:,2) = idvg_p(:,2);
idvg_p(:,3) = idvg_p(:,3);

end
