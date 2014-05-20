function [timeCourse,population] = stochastic(duration,cI_r0,cI_p0,cro_r0,cro_p0)
	cI_degrade = 0.8; % Equivalent for RNA and protein
	cro_degrade = 1.2; % Equivalent for RNA and protein
	formation = 50; % Equivalent for cro and cI; both proteins and RNA
	halfValue = 10; % Equivalent for cro and cI

	% Create a vector to store the time course
    timeCourse = zeros(1,duration);
    
    %{
    Create a vector to store the populations
               1_2 3_4 5_6 7_8 9_... n
    cI_rna    |
    cI_prot   |
    cro_rna   |
    cro_prot  |
    %}
   population = [timeCourse;timeCourse;timeCourse;timeCourse];
   population(1,1) = cI_r0;
   population(2,1) = cI_p0;
   population(3,1) = cro_r0;
   population(4,1) = cro_p0;
   
   % Pre-allocate a series of random integers for us to use. This is
   % faster.
   random = rand(2,duration);
   
   % Define the state changes depending on the reaction that occurs
   % 4 rows for each species. 8 columns. Each column refers to the reaction.
   stateChanges = [
       0 0 1 -1 0 0 0 0; % cI_rna
       1 -1 0 0 0 0 0 0; % cI_prot
       0 0 0 0 0 0 1 -1; % cro_rna
       0 0 0 0 1 -1 0 0  % cro_prot
   ];
   
   for i=1:duration
       propensity = zeros(1,8);
       propensity(1) = formation * population(1,i);
       propensity(2) = cI_degrade * population(2,i);
       propensity(3) = formation*(1-((population(4,i)^2)/(halfValue^2 + population(4,i)^2)));
       propensity(4) = cI_degrade * population(1,i);
       propensity(5) = formation * population(3,i);
       propensity(6) = cro_degrade * population(4,i);
       propensity(7) = formation*(1-((population(2,i)^2)/(halfValue^2 + population(2,i)^2)));
       propensity(8) = cro_degrade * population(3,i);
       
       a0 = sum(propensity);
       
       deltaT = -log(random(1,i))./a0;
       reactionIndex = find(cumsum(propensity) >= random(2,i)*a0 , 1 , 'first');
       
       population(:,i+1) = max(population(:,i) + stateChanges(:,reactionIndex),0);
       timeCourse(i+1) = timeCourse(i) + deltaT;
   end
end