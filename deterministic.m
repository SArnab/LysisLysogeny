function [cI_r,cI_p,cro_r,cro_p] = deterministic(duration,timestep,cI_r0,cI_p0,cro_r0,cro_p0)
	cI_degrade = 1.2; % Equivalent for RNA and protein
	cro_degrade = 0.8; % Equivalent for RNA and protein
	formation = 50; % Equivalent for cro and cI; both proteins and RNA
	halfValue = 10; % Equivalent for cro and cI

	% Create time-dependent vectors to store values
	timeVector = (0:timestep:duration);
	cI_r = zeros(1,length(timeVector));
	cI_p = zeros(1,length(timeVector));
	cro_r = zeros(1,length(timeVector));
	cro_p = zeros(1,length(timeVector));
	% Declare initial concentrations
	cI_r(1) = cI_r0;
	cI_p(1) = cI_p0;
	cro_r(1) = cro_r0;
	cro_p(1) = cro_p0;

	for i=1:length(timeVector)-1
		cI_p(i+1) = cI_p(i) + timestep*(formation*cI_r(i) - cI_degrade*cI_p(i));
		cI_r(i+1) = cI_r(i) + timestep*(formation*(1-(cro_p(i)^2 / (halfValue^2 + cro_p(i)^2))) - cI_degrade*cI_r(i));

		cro_p(i+1) = cro_p(i) + timestep*(formation*cro_r(i) - cro_degrade*cro_p(i));
		cro_r(i+1) = cro_r(i) + timestep*(formation*(1-(cI_p(i)^2 / (halfValue^2 + cI_p(i)^2))) - cro_degrade*cro_r(i));
	end

	plot(timeVector,cI_r);
end