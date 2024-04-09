%% COM DSA MODEL version 4 used in the BLOG%%
function runDsaStochModel4_blog(sfa_method,scenario,apply_debt_safeguard,...
    plotting,power,plausibility,language,stoch_method,saveFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% SELECT SCENARIO AND SFA METHOD %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     SELECT SFA METHOD:
%          0) COM ZERO ASSUMPTION
%         -1) ALTERNATIVE ASSUMPTION      
%     SELECT SCENARIO:
%         1) ADJUSTMENT
%         2) LOWER SPB
%         3) ADVERSE r-g
%         4) FINANCIAL STRESS
%     APPLY DEBT SAFEGUARD
%         yes=1, no=0
%     PLOTTING
%         yes=1, no=0
%     Stochastic samples (power)
%         1 = 10 simulated paths
%         2 = 100
%         3 = 1000
%         4 = 10 0000
%         5 = 100 0000
%         6 = million
%     Plausibility value
%         7 = 70%
%         8 = 80%
%         9 = 90%
%     Language for Stochastic plots
%         1 = English
%         2 = Suomi
%         3 = Swedish
%     Stochastic Method
%         1 = Normal
%         2 = Bootstrap
%     Save results as -mat (OPTIONAL)
%         1 = save
% Code calculates minimum yearly adjustment for four year adjustment plan
% for Finland
%
% Code produces debt projections following the structure of 
% the European Commission's Debt Sustainability Monitor 2023. 
%
% Currently, all four scenarios are available: baseline adjustment
% scenario, lower structural primary balance (SPB) scenario, 
% Adverse r-g scenario and Financial stress sccenario. Also, 
% stock-flow-adjustment (SFA) method can be chosen to follow 
% Commission assumption, or an alternative assumption based on 
% linearly declining SFA.
%
% TIMING: First observation corresponds to year 2022 in data matrix
%           the years run from 2022 to 2038 (17 obs). 2022-2024 are 
%           pre-adjustment plan periods, 2025-2028 are adjustment plan
%           periods and 2029-2038 are post-adjustment plan periods.
%
% NOTE: Code uses function project_debt4v.m to find minimum 
%           yearly adjustment which satisfy only criteria A
%
% For comments and suggestions please contact peetu.keskinen[at]vtv[dot]fi
% 
% Author: Peetu Keskinen
% Date: 27/2/2024
%
%% MAKE NOTATION MORE GENERAL!


%% DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step_size=0.01; % define step size for adjustment
a = 0:step_size:1;  % Grid of possible values of a
deter = 1; % use this to project debt using deterministic scenarios

% constants related to the adjustment phase 
adj_periods = 4; % number of adjustment periods (2025 - 2028)
pre_plan_periods = 3;   % number of periods before the adjustment plan (22-24)
adjustment_start = pre_plan_periods + 1; % start period of adjustment (2025)
adjustment_end = pre_plan_periods + adj_periods; % end period of adjustment
post_plan_periods = 10; % number of periods after the adjustment plan (2028-38)
total_periods = pre_plan_periods + adj_periods + post_plan_periods;

fiscal_periods = 2; % number of periods in which adjustment effects 
                          %     output gap beyond adjustment plan
adjusted_post_plan_periods = post_plan_periods - fiscal_periods;

ameco_end_t = 3; %the last period of AMECO forecast (2024)
sfa_periods = 10; % number of periods to decrease sfa during projection
%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data from the file
% Important to choose right range!
data = readmatrix('AmecoFinlandData4v_1.xlsx', 'Range', 'B49:V65','Sheet','COM');
BBGdata = readmatrix('AmecoFinlandData4v_1.xlsx', 'Range', 'K4:N20','Sheet','BBG');
ECBdata = readmatrix('AmecoFinlandData4v_1.xlsx', 'Range', 'B3:G3','Sheet','ECB');

% Name variables 
%iir = data(3:end,14)./100; % implicit interest rate
potgdp = data(1:end,6); % potential gdp level
inflation = data(1:end,10)./100; % inflation rate
og = data(1:end,8); % output gap
dcoa=data(1:end,19); % delta cost of ageing reference year 2028

debt_ltn = zeros(total_periods, 1);
debt_total = zeros(total_periods, 1);
debt_total(1) = ECBdata(2); % total debt
debt_st = zeros(total_periods, 1);
debt_st(1) = ECBdata(3); % short-term debt
debt_lt = zeros(total_periods, 1);
debt_lt(1) = ECBdata(4); % long-term debt
share_lt_maturing_t0 = ECBdata(5); % current share of lt maturing debt (2022)
share_lt_maturing_t10 = ECBdata(6); % target share of lt maturing debt (2032)
% current value approaches target value
share_lt_maturing = linspace(share_lt_maturing_t0,share_lt_maturing_t10,11)';
share_lt_maturing = [share_lt_maturing; share_lt_maturing_t10*ones(6,1)];
iir = zeros(total_periods,1);
iir(1:3) = BBGdata(1:3,4); % use 2022-2024 COM forecasts for IIR
i_st = BBGdata(1:end,2); % short-term market interest rate
i_lt = BBGdata(1:end,3); % long-term market interest rate

%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Commission fiscal multiplier based on Carnot and de Castro (2015)
phi = 0.75; % COM assumption on fiscal multiplier

% Fiscal policy effects output gap during the consolidation 
% period 2025 - 2028 + 2 years (see 2nd infobox p. 66 from the report)
  
if adj_periods == 4 % 4 year plan (2 peak multipliers)
    m = [zeros(1,pre_plan_periods) phi phi*(5/3) phi*(6/3) phi*(6/3)...
      phi phi*(1/3) zeros(1,adjusted_post_plan_periods)]';
  
elseif adj_periods == 7 % 7 year plan (5 peak multipliers)
    m = [zeros(1,pre_plan_periods) phi phi*(5/3) phi*(6/3) phi*(6/3)...
        phi*(6/3) phi*(6/3) phi*(6/3) phi phi*(1/3) ...
        zeros(1,adjusted_post_plan_periods)]';
end

epsilon = 0.582; %semi-elasticity of budget balance for Finland
debt_initial = 73.2855; %debt in 2022, source:ECB
%debt_initial = 74.3; %debt in 2023, COM forecast 2023 autumn

spb = NaN(total_periods,1);
spb(1:3) = data(1:3,12);    % COM estimate of spb 22-24
pb = NaN(total_periods,1);
pb(1:3) = data(1:3,13);     % COM estimate of pb 22-24
rgdp = NaN(total_periods,1);
rgdp(1:3) = data(1:3,4);    % COM estimate of real gdp 22-24

PB = zeros(total_periods, 1); % Primary balance in euros
PB(1:3) = data(1:3,20);       % COM estimate
SF = zeros(total_periods, 1); % Stock flow adjustment in euros
SF(1:3) = data(1:3,21);       % COM estimate

% Create Debt Sustainability Safequard
if debt_initial>=60 && debt_initial<=90
    debt_safeguard=-0.5; % debt must decline 0.5 %-points per year
elseif debt_initial>90
    debt_safeguard=-1; % debt must decline 1.0 %-points per year
else
    disp('Debt below 60%, Debt Sustainability safeguard not needed.');
end

%% SFA METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COM assumption
if sfa_method == 0    % stock-flow adjustment with t+2 on sfa=0
    
        sfa=data(1:end,16); % sfa data 2022-2038

% Alternative assumption
elseif sfa_method == -1       % linear sfa correction
        sfa=data(1:end,16);   % sfa data 2022-2038

for j= ameco_end_t : ameco_end_t + sfa_periods - 1
    
        % decrease sfa linearly
        sfa(j+1)=sfa(j) - sfa(ameco_end_t)/sfa_periods; 
end
end

%% CALCULATE DEBT PATHS FOR ALL VALUES OF ADJUSTMENT

% display selected scenario and if safeguard apply
disp(['Selected scenario is ', num2str(scenario)]);
disp(['Does Safeguard apply? Yes: 1, No: 0. Selection is ', num2str(apply_debt_safeguard)]);

% Loop over each a value of adjustment in the grid
D = zeros(total_periods,length(a)); % shell for debt
Gn = zeros(total_periods,length(a)); % shell for nominal growth
Gr = zeros(total_periods,length(a)); % shell for real growth

% Calculate paths for ALL values of a
for i = 1:length(a)
    [D(:,i),Gn(:,i),Gr(:,i),~,~,~] = project_debt4v(scenario,a(i), iir, potgdp,...
        og, epsilon,m,dcoa,sfa,inflation, debt_initial, spb,...
        rgdp,debt_st,debt_total,debt_ltn,debt_lt,i_st,i_lt,...
        share_lt_maturing,PB,SF,pb,deter,zeros(total_periods,1),...
        zeros(total_periods,1),zeros(total_periods,1));
end

%Start checking the minimum adjustment until declining=true
    
solution_column = 0;  % Initialize to 0 to indicate no solution found

%% CHECK DECLINING DEBT AND POSSIBLE SAFEGUARD
if apply_debt_safeguard == 1
% Declining debt path in 10 year review period and Debt sustainability safeguard
for j = 1:length(a)  % Loop through columns of D
    if all(diff(D(adjustment_end:end,j)) < 0) && ... 
          mean(diff(D(adjustment_start-1:adjustment_end,j))) < debt_safeguard 
        
        solution_column = j;
        break;  % Exit the loop once the condition is met
    end
end

elseif apply_debt_safeguard == 0
    % Declining debt path in 10 year review period
for j = 1:size(D,2)  % Loop through columns of D
    if all(diff(D(adjustment_end:end,j)) < 0)
        
        solution_column = j;
        break;  % Exit the loop once the condition is met
    end
end
end

if solution_column == 0
    disp('No value could be found that satisfies the condition.');
else
    optimal_a = a(solution_column); % minimum consolidation a*
    disp(['Minimum required adjustment (Deterministic) is a*=', num2str(optimal_a)]);
end

%% PROJECT DEBT USING MINIMUM ADJUSTMENT a*
    [debt_path,optimal_Gn,optimal_Gr,iir_path,pb_path,spb_path] = project_debt4v(scenario,optimal_a, iir, potgdp,...
        og, epsilon,m,dcoa,sfa,inflation, debt_initial, spb,...
        rgdp,debt_st,debt_total,debt_ltn,debt_lt,i_st,i_lt,...
        share_lt_maturing,PB,SF,pb,deter,zeros(total_periods,1),...
        zeros(total_periods,1),zeros(total_periods,1));
    
disp(['Structural primary balance in 2029 is SPB*=', num2str(round(spb_path(adjustment_end),2))]);   
disp(['Structural primary balance in 2039 is SPB*=', num2str(round(spb_path(end),2))]);   

%% START PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotting == 1
    %% Plot 2D the debt path and related variables
    FigSize = 600; %set figure size
    time =(2022:2038)';
    figure(1);
    set(gca,'fontname','Calibri') 
    subplot(1,3,1);
    plot(time, debt_path, '-o', 'LineWidth', 2,'Color','#002C5F');
    title('(a)','FontSize',20,'FontName', 'Calibri')
    ylim([60 max(debt_path)+10])
    ylabel('Velkasuhde, %','FontSize',20,'FontName', 'Calibri');
    legend('Velka/BKT','Location','southeast','FontName', 'Calibri'...
        ,'FontSize',8);
    grid on;

    subplot(1,3,2);
    plot(time, spb_path, '-o', 'LineWidth', 2,'Color','#002C5F');
    hold on
    plot(time, pb_path, '--', 'LineWidth', 2,'Color','#FF6875');
    title('(b)','FontSize',20,'FontName', 'Calibri')
    ylabel('% suhteessa BKT:hen','FontSize',20,'FontName', 'Calibri');
    legend({"Rakenteellinen perusjäämä" + newline + "(ikäsidonnaisilla menoilla korjattu)",'Perusjäämä'},'Location','southeast',...
        'FontSize',8,'FontName', 'Calibri');
    grid on;

    subplot(1,3,3);
    plot(time, 100*optimal_Gn, '-o', 'LineWidth', 2,'Color','#002C5F');
    hold on
    plot(time, iir_path, '--', 'LineWidth', 2,'Color','#FF6875');
    plot(time, 100*optimal_Gr, 'LineWidth', 2,'Color','#FDB74A');
    title('(c)','FontSize',20,'FontName', 'Calibri')
    %title('$g_{t}$ \& $r_{t}$','interpreter','latex',...
    %   'FontSize',20,'FontName', 'Calibri')
    ylabel('vuosikasvu, %','FontSize',20,'FontName', 'Calibri');
    legend('Nimellinen BKT-kasvu','Nimellinen Korkotaso','Reaalinen BKT-kasvu','Location',...
        'southeast','FontName', 'Calibri','FontSize',8);
    grid on;
    if scenario==1
    sgtitle('Perusura','FontSize',24,'FontName', 'Calibri');
    elseif scenario==2
    sgtitle('Epäsuotuisa SPB','FontSize',24,'FontName', 'Calibri');
    elseif scenario==3
    sgtitle('Epäsuotuisa r-g','FontSize',24,'FontName', 'Calibri');
    elseif scenario==4
    sgtitle('Rahoitusmarkkinahäiriö','FontSize',24,'FontName', 'Calibri');
    end
    % Set the figure size
    set(gcf, 'Position', [20 20 2*FigSize FigSize]);

    % Save figure data %
    if scenario==1
    PlotDataAdj=[time debt_path spb_path pb_path...
                100*optimal_Gn 100*iir_path 100*optimal_Gr];
    % Assuming PlotDataLowerAdj is a matrix with the correct number of columns to match the headers.
    headers = {'vuosi', 'velkasuhde', 'rakenteellinen perusjaama',...
        'perusjaama', 'nimellinen bkt kasvu', 'nimellinen korkotaso', 'reaalinen bkt kasvu'};
    % Convert the matrix to a table
    T = array2table(PlotDataAdj, 'VariableNames', headers);
    % Write the table to a text file with headers
    writetable(T, "PlotDataAdj.txt", 'Delimiter', ' ');
    print('Adjustment','-dpng', '-r300','-cmyk');

    elseif scenario==2
    PlotDataLowerSPB=[time debt_path spb_path pb_path...
                100*optimal_Gn 100*iir_path 100*optimal_Gr];
    % Assuming PlotDataLowerSPB is a matrix with the correct number of columns to match the headers.
    headers = {'vuosi', 'velkasuhde', 'rakenteellinen perusjaama',...
        'perusjaama', 'nimellinen bkt kasvu', 'nimellinen korkotaso', 'reaalinen bkt kasvu'};
    % Convert the matrix to a table
    T = array2table(PlotDataLowerSPB, 'VariableNames', headers);
    % Write the table to a text file with headers
    writetable(T, "PlotDataLowerSPB.txt", 'Delimiter', ' ');
    print('LowerSPB', '-dpng', '-r300','-cmyk');
    
        elseif scenario==3
    PlotDataAdverseR_G=[time debt_path spb_path pb_path...
                100*optimal_Gn 100*iir_path 100*optimal_Gr];
    % Assuming PlotDataArverseR_G is a matrix with the correct number of columns to match the headers.
    headers = {'vuosi', 'velkasuhde', 'rakenteellinen perusjaama',...
        'perusjaama', 'nimellinen bkt kasvu', 'nimellinen korkotaso', 'reaalinen bkt kasvu'};
    % Convert the matrix to a table
    T = array2table(PlotDataAdverseR_G, 'VariableNames', headers);
    % Write the table to a text file with headers
    writetable(T, "PlotDataAdverseR_G.txt", 'Delimiter', ' ');
    print('AdverseR_G', '-dpng', '-r300','-cmyk');
    
        elseif scenario==4
    PlotDataFinancialStress=[time debt_path spb_path pb_path...
                100*optimal_Gn 100*iir_path 100*optimal_Gr];
    % Assuming PlotDataFinancialStress is a matrix with the correct number of columns to match the headers.
    headers = {'vuosi', 'velkasuhde', 'rakenteellinen perusjaama',...
        'perusjaama', 'nimellinen bkt kasvu', 'nimellinen korkotaso', 'reaalinen bkt kasvu'};
    % Convert the matrix to a table
    T = array2table(PlotDataFinancialStress, 'VariableNames', headers);
    % Write the table to a text file with headers
    writetable(T, "PlotDataFinancialStress.txt", 'Delimiter', ' ');
    print('FinancialStress', '-dpng', '-r300','-cmyk');
    end

    %% Plot 3D plots (real gdp growth)
    % Define the colors in RGB, normalized to [0, 1]
    color1 = [204, 213, 223] / 255; % Light blue 2
    color2 = [153, 171, 191] / 255; % Light blue 1
    color3 = [102, 128, 159] / 255; % Medium blue 2
    color4 = [51, 86, 127] / 255;   % Medium blue 1
    color5 = [0, 44, 95] / 255;     % Dark blue

    % Preallocate an array for the colormap
    numColors = 10;
    customColormap = zeros(numColors, 3);

    % Generate intermediate colors using linspace
    for i = 1:3 % For each color channel
        % Combine all five colors' current channel into an array
        originalChannels = [color1(i), color2(i), color3(i), color4(i), color5(i)];
        % Interpolate to find two intermediate colors between each original color
        customColormap(:, i) = interp1(1:length(originalChannels), originalChannels, linspace(1, length(originalChannels), numColors));
    end

    % Apply the custom colormap to the current figure
    colormap(customColormap);
    figure(2);
    %subplot(1,2,1);
    ax1 = gca;  % Get handle to current axes
    x1 = time;
    y1 = a;
    [X1, Y1] = meshgrid(x1, y1);
    Z1 = 100 * Gr;
    surf(Y1', X1', Z1, 'FaceAlpha', 1);
    xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 20)
    zlabel('BKT-kasvu, %', 'FontSize', 20, 'FontName', 'Calibri');
    zlim([0 max(100*Gr,[], 'all') + 1])
    colormap(ax1, customColormap); 
    set(gcf, 'Position', [250 250 FigSize FigSize]);

    %% Plot 3D plots (debt-to-gdp ratio)
    figure(3)
    ax2 = gca;  % Get handle to current axes
    x2 = time;
    y2 = a;
    [X2, Y2] = meshgrid(x2, y2);
    Z2 = D;
    surf(Y2', X2', Z2, 'FaceAlpha', 1);
    xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 20)
    zlabel('Velka/BKT', 'FontSize', 20, 'FontName', 'Calibri');
    view( 156.5618,  12.3745); %adjust view angle
    zlim([0  max(D(end,:))+10])
    colormap(ax2, customColormap);  % Apply custom colormap to second subplot

    % Set color data based on Z-values (debt ratio change)
    caxis(ax2, [40, 90]);
    % Add colorbar
    colorbar('Position', [0.95, 0.17, 0.02, 0.8], 'FontSize', 14, 'FontName', 'Calibri');

    % Set the font name for all text objects in the current figure
    set(gca, 'FontName', 'Calibri');       % Change font for axes tick labels
    set(findall(gcf,'type','text'), 'FontName', 'Calibri'); % Change font for titles, labels, legends, etc.

    % Change the font size as well
    set(gca, 'FontSize', 20);            % Change font size for axes tick labels
    set(findall(gcf,'type','text'), 'FontSize', 20); % Change font size for titles, labels, legends, etc.

    % Set the figure size
    set(gcf, 'Position', [100 500 1.5*FigSize 1.5*FigSize]);

    if scenario==1
    % Increase the resolution to 300 dpi
    print('Adjustment3D','-dpng', '-r300','-cmyk');
    elseif scenario==2
    % Increase the resolution to 300 dpi
    print('Lower3D','-dpng', '-r300','-cmyk');
    elseif scenario==3
    % Increase the resolution to 300 dpi
    print('AdverseR_G3D','-dpng', '-r300','-cmyk');
    elseif scenario==4
    % Increase the resolution to 300 dpi
    print('FinancialStress3D','-dpng', '-r300','-cmyk');
    end

    %% Bird view plot (debt-to-gdp ratio)
    figure(4)
    ax2 = gca;  % Get handle to current axes
    x2 = time;
    y2 = a;
    [X2, Y2] = meshgrid(x2, y2);
    Z2 = D;
    surf(Y2', X2', Z2, 'FaceAlpha', 1);
    xlabel('$a$', 'Interpreter', 'latex', 'FontSize', 20)
    zlabel('Velka/BKT', 'FontSize', 20, 'FontName', 'Calibri');
    view(2); % Adjust th view angle
    zlim([0  max(D(15,:))+10])
    colormap(ax2, customColormap);  % Apply custom colormap to second subplot

    % Set color data based on Z-values (debt ratio change)
    caxis(ax2, [40, 90]);
    colorbar

    if scenario==1
    sgtitle('Perusura','FontSize',30,'FontName', 'Calibri');
    elseif scenario==2
    sgtitle('Epäsuotuisa SPB','FontSize',30,'FontName', 'Calibri');
    elseif scenario==3
    sgtitle('Epäsuotuisa r-g','FontSize',30,'FontName', 'Calibri');
    elseif scenario==4
    sgtitle('Rahoitusmarkkinahäiriö','FontSize',30,'FontName', 'Calibri');
    end

    % Set the figure size
    set(gcf, 'Position', [200 600 FigSize FigSize]);

    if scenario==1
    % Increase the resolution to 300 dpi
    print('AdjustmentBird','-dpng', '-r300','-cmyk');
    elseif scenario==2
    % Increase the resolution to 300 dpi
    print('LowerBird','-dpng', '-r300','-cmyk');
    elseif scenario==3
    % Increase the resolution to 300 dpi
    print('AdverseR_GBird','-dpng', '-r300','-cmyk');
    elseif scenario==4
    % Increase the resolution to 300 dpi
    print('FinancialStressBird','-dpng', '-r300','-cmyk');
    end

elseif plotting==0
    disp(['*Plotting not selected*']);
else
    disp(['Define plotting variable']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% STOCHASTIC SCENARIO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load quarterly first differenced data 1999Q2-2022Q4
% data in terms of %-points of GDP
% Columns: short-term interest rate (%-points), 
%          long-term interest rate (%-points), 
%          nominal gdp growth (%-points), 
%          primary balance (%-points)

dataStoch = readmatrix('AmecoFinlandData4v_1.xlsx', 'Range', 'I3:L97','Sheet','STOCH');
dataStochBoot = readmatrix('AmecoFinlandData4v_1.xlsx', 'Range', 'C3:F49','Sheet','STOCH'); 
% C3 (1976)/ C17 (1990) / C26 1999
alpha = debt_st(1) / debt_total(1); % share of short-term debt
meanValues = mean(dataStoch);
stdDevValues = std(dataStoch);
% correct for outliers
outlier_threshold = 3;
lowerThreshold = meanValues - outlier_threshold * stdDevValues;
upperThreshold = meanValues + outlier_threshold * stdDevValues;
% Clip values based on the thresholds
dataClipped = max(dataStoch, lowerThreshold);
dataClipped = min(dataClipped, upperThreshold);
% Create a logical mask for replaced values
replacedMask = (dataClipped ~= dataStoch);
% Calculate the number of replaced values
numReplacedValues = nnz(replacedMask);
disp(['Number of replaced values: ', num2str(numReplacedValues)]);

%% Parameters
share_lt_maturing_t10 = 0.039528;% ECB data, REMOVE ONCE MERGED
nbr_q = 4; % number of quarters in a year
T_stochastic = 5; % years affected by stochastic shocks
nbr_q_shocks = 20*10^power; % number of quarterly shocks to generate
post_stoch = 5; % post stochastic years in evaluation phase
pre_stoch = 7; % pre stochastic years before the plan

%% Definitions
nbr_vars = size(dataClipped,2); % number of variables
nbr_q_gen = nbr_q*T_stochastic; % quarters affected by stochastic shocks
nbr_sim_paths = nbr_q_shocks/nbr_q_gen; % number of 5-year forecasts (sim. debt paths)
nbr_y_shocks = nbr_q_shocks/nbr_q; % number of years in the forecast 
% Average maturity
m_res_lt = 1/share_lt_maturing_t10; % average residual maturity of lt bonds
maturity_quarters = nbr_q*m_res_lt; % average residual maturity in quarters

% Name methods
if stoch_method == 1
    methodName = 'normal';
elseif stoch_method == 2
    methodName = 'boot';
else
    methodName = 'unknown';
end

%% Generate quarterly shocks
% Sampling methods
if stoch_method == 1
    % Normal Random Sampling
    mu = zeros(1, nbr_vars);
    sigma = cov(dataClipped);
    randSamples = mvnrnd(mu, sigma, nbr_q_shocks);
    Sample = randSamples; 

    % reshape into yearly 4x4 matrices. now row: vars, column: quarter)
    Rands=reshape(Sample',[nbr_q,nbr_vars,nbr_y_shocks]); % reshape

    %% Aggregate quarterly shocks to yearly
    % create shell for yearly shocks
    Shock = zeros(nbr_y_shocks,nbr_vars); 
    % store yearly shocks
    e_i_st = zeros(nbr_sim_paths,T_stochastic); % shells
    e_g = zeros(nbr_sim_paths,T_stochastic);
    e_pb = zeros(nbr_sim_paths,T_stochastic);

    % get yearly shock by summing quarterly shocks
    for i=1:nbr_y_shocks
    Shock(i,:)=sum(Rands(:,:,i),2)';
    end

    % reshape (variable, 5-year forecast, N copies)
    Shock = reshape(Shock',[nbr_vars T_stochastic nbr_sim_paths]);

    for i=1:nbr_sim_paths
    e_i_st(i,:) = Shock(1,:,i); % short term market interest rate
    e_g(i,:) = Shock(3,:,i);    % nominal gdp growth rate
    e_pb(i,:) = Shock(4,:,i);   % primary balance
    end

%% Construct long term interest rate shocks
    shock_i_lt = Sample(:,2); % generated quarterly i_lt shocks
    % use function sumq2y to get persistent yearly shocks
    [e_i_lt] = sumq2y(shock_i_lt,m_res_lt,T_stochastic);

    % add zeros to get correct dimensions
    % iir shock
    i_st_shock = [zeros(pre_stoch,nbr_sim_paths); e_i_st'; zeros(post_stoch,nbr_sim_paths)];
    i_lt_shock = [zeros(pre_stoch,nbr_sim_paths); e_i_lt'; zeros(post_stoch,nbr_sim_paths)];

    iir_shock = alpha.*i_st_shock + (1-alpha).*i_lt_shock; 
    % gdp and pb shocks
    g_shock = [zeros(pre_stoch,nbr_sim_paths); e_g'; zeros(post_stoch,nbr_sim_paths)];
    pb_shock = [zeros(pre_stoch,nbr_sim_paths); e_pb'; zeros(post_stoch,nbr_sim_paths)];

%% Project debt paths under stochastic scenarios

    % shell for final 3D stochastic debt matrix 
    D_stoch = zeros(total_periods,length(a),nbr_sim_paths);

    % Initialize waitbar
    h = waitbar(0, 'Running simulations...');

    % Loop through stochastic shocks...
    for k = 1:nbr_sim_paths

        % Update waitbar
        waitbar(k/nbr_sim_paths, h, sprintf('Calculate stochastic debt path %d of %d', k, nbr_sim_paths));

        D_temp = zeros(total_periods,length(a)); % shell
        iir_stoch = iir_shock(:,k); % iir
        g_stoch = g_shock(:,k); % nominal gdp growth
        pb_stoch = pb_shock(:,k); % primary balance

        % ... and through values of a
        for i = 1:length(a)

            [D_temp(:,i),~,~,~,~,~] = project_debt4v(scenario,a(i), iir, potgdp,...
                og, epsilon,m,dcoa,sfa,inflation, debt_initial, spb,...
                rgdp,debt_st,debt_total,debt_ltn,debt_lt,i_st,i_lt,...
                share_lt_maturing,PB,SF,pb,stoch_method,g_stoch,pb_stoch,iir_stoch);
        end

    % store 2D matrix
    D_stoch(:,:,k) = D_temp;
    end

    % Close waitbar after completion
    close(h);

elseif stoch_method == 2
    
    % shell for final 3D stochastic debt matrix 
    D_stoch = zeros(total_periods,length(a),nbr_sim_paths);
    % Sampling with Block-Bootstrap Approach
    %blockSize = 2; % Size of each block (2 years2)
    nbr_boots = 5; % number of bootstrap draws in one simulation
    
    % load Stoch Boot data (d_rgdp	pb	iir	sfa)
    sampleSize = size(dataStochBoot, 1); % sample size
    % Initialize bootstrapSamples as a 3D array
    bootstrapSamples = zeros(total_periods, size(dataStochBoot, 2), nbr_sim_paths);
    blockStartIdx = zeros(nbr_boots, 1); % Initialize the vector
    
    % create block bootstrap sample
    for i = 1:nbr_sim_paths
    
        for j = 1:2:nbr_boots
            % Draw a random index for odd elements
            blockStartIdx(j) = randi([1, sampleSize - 1]);

            % For even elements, add 1 to the previous odd element
            if j+1 <= nbr_boots 
                blockStartIdx(j+1) = blockStartIdx(j) + 1;
            end
        end

        % Select the blockSize number of consecutive years
        selectedBlock = dataStochBoot(blockStartIdx, :);
        bootStart = adjustment_end+1;
        bootEnd = adjustment_end+nbr_boots;
        
        % Store is in 3D zero metrix
        bootstrapSamples(bootStart:bootEnd,:,i) = selectedBlock;
        
    end
   
    
    % Initialize waitbar
    h = waitbar(0, 'Running simulations...');
    
    
        % Loop through stochastic shocks...
        for k = 1:nbr_sim_paths
        
        % Update waitbar
        waitbar(k/nbr_sim_paths, h, sprintf('Calculate stochastic debt path %d of %d', k, nbr_sim_paths));

        D_temp = zeros(total_periods,length(a)); % shell   
        % pick variables  (d_rgdp	pb	iir	sfa)
        g_stoch = bootstrapSamples(:,1,k) ; % nominal gdp growth
        pb_stoch = bootstrapSamples(:,2,k) ; % primary balance
        iir_stoch = bootstrapSamples(:,3,k) ; % iir
        %sfa_stoch = bootstrapSamples(:,4,k) ; % sfa
        
        % ... and through values of a
        for i = 1:length(a)

            [D_temp(:,i),~,~,~,~,~] = project_debt4v(scenario,a(i), iir, potgdp,...
                og, epsilon,m,dcoa,sfa,inflation, debt_initial, spb,...
                rgdp,debt_st,debt_total,debt_ltn,debt_lt,i_st,i_lt,...
                share_lt_maturing,PB,SF,pb,stoch_method,g_stoch,pb_stoch,iir_stoch);
        end

        % store 2D matrix in a 3D matrix
        D_stoch(:,:,k) = D_temp;
        end

    % Close waitbar after completion
    close(h);    

    else
        error('Invalid selection. Please choose 1 or 2.');
    
end


%% STOCHASTIC DEBT PATHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% Initialize variables and load data
PlanPeriods=4;
prePlanPerids=3;
horizon=T_stochastic+PlanPeriods+prePlanPerids;

figure;
hold on;
FigSize = 600;
range = 10:10:90; % percentile range 10% to 90%
% Empty shell (horizon, number of percentiles, number of adjustment values)
debt_path_stoch = zeros(horizon, length(range), length(a));

% calculate percentiles
for k = 1:length(a)
    for j = 1:12
        debt_path_stoch(j, :,k) = prctile(D_stoch(j,k, :), range);
    end
end

%Start checking the minimum adjustment until plausibility ok, then declining=true
solution_columnZ = 0;  % Initialize to 0 to indicate no solution found

%% CHECK DECLINING DEBT AND POSSIBLE SAFEGUARD
if apply_debt_safeguard == 1
% Check the plausibility criteria
for m = 1:length(a)  % Loop through values of a
    % 1) must be less in the last period than after the adjustment period
    % AND 2) debt safegueard must be satisfied
    
    if debt_path_stoch(horizon,plausibility,m) < ...
            debt_path_stoch(adjustment_end,plausibility,m) &&...
        mean(diff(debt_path_stoch(adjustment_start-1:adjustment_end,plausibility,m)))...
        < debt_safeguard 
        
        solution_columnZ = m;
        debt_path_success = debt_path_stoch(:,:,solution_columnZ);
        break;  % Exit the loop once the condition is met
    end
end

elseif apply_debt_safeguard == 0
% Check the plausibility criteria
for l = 1:length(a)  % Loop through values of a
    % must be less in the last period than after the adjustment period
    if debt_path_stoch(horizon,plausibility,l) < debt_path_stoch(adjustment_end,plausibility,l)
        
        solution_columnZ = l;
        debt_path_success = debt_path_stoch(:,:,solution_columnZ);
        break;  % Exit the loop once the condition is met
    end
end
end

if solution_columnZ == 0
    disp('No value could be found that satisfies the condition.');
else
    optimal_stoch_a = a(solution_columnZ); % minimum consolidation a*
    disp(['Minimum required adjustment (Stochastic) is a*=', num2str(optimal_stoch_a)]);
end

%% Plot the results
% Change the axes' properties
set(gca, 'FontName', 'Calibri','FontSize',12)

% Plot the results
time_period = 2022:2033;

% % Define the function to fill areas between percentiles
% fillArea = @(p1, p2, color) fill([time_period, fliplr(time_period)],...
%              [debt_path_success(:,p1)', fliplr(debt_path_success(:,p2)')],...
%                 color, 'LineStyle', 'none');

            % NEW STARTS
% Modify the fillArea function to include an invisible marker plot
fillArea = @(p1, p2, color) ...
    plotFillAndMarker(time_period, debt_path_success(:, p1)', debt_path_success(:, p2)', color);

% New function to plot fill and invisible marker
function h = plotFillAndMarker(x, y1, y2, color)
    fill([x, fliplr(x)], [y1, fliplr(y2)], color, 'LineStyle', 'none'); % Fills the area between curves    
    % Plots visible markers outside the plot area for legend purposes
    h = plot(nan, nan, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
end


            % NEW ENDS



% Fill areas for different percentile ranges
% % VTV colors
h1 = fillArea(1, 9, [153/255 171/255 191/255]); % 10-90th percentile, lightest blue
h2 = fillArea(2, 8, [102/255 128/255 159/255]); % 20-80th percentile, intermediate blue
h3 = fillArea(3, 7, [51/255 86/255 127/255]); % 30-70th percentile, darker blue 
h4 = fillArea(4, 6, [0/255 44/255 95/255]); % 40-60th percentile, darkest blue 

% VTV colors
% fillArea(1, 9, [204/255 213/255 223/255]); % 10-90th percentile, lightest blue '#CCD5DF'
% fillArea(2, 8, [153/255 171/255 191/255]); % 20-80th percentile, intermediate blue '#99ABBF'
% fillArea(3, 7, [102/255 128/255 159/255]); % 30-70th percentile, darker blue '#66809F'
% fillArea(4, 6, [51/255 86/255 127/255]); % 40-60th percentile, darkest blue '#33567F'
% % Gray colors
% fillArea(1, 9, [0.9 0.9 0.9]); % 10-90th percentile, lightest grey
% fillArea(2, 8, [0.7 0.7 0.7]); % 20-80th percentile, intermediate grey
% fillArea(3, 7, [0.5 0.5 0.5]); % 30-70th percentile, darker grey
% fillArea(4, 6, [0.3 0.3 0.3]); % 40-60th percentile, darkest grey

% Add a red horizontal line at the year 2024 and 2028 
yLimits = get(gca, 'ylim'); % Get the current y-axis limits
% VTV red
h5 = line([2024 2024], yLimits, 'Color','#fdb84a', 'LineStyle', '-', 'LineWidth', 1);
h6 = line([2028 2028], yLimits, 'Color','#fdb84a', 'LineStyle', '--', 'LineWidth', 1.5);
% regular red
% line([2024 2024], yLimits, 'Color', 'red', 'LineStyle', '-');
% line([2028 2028], yLimits, 'Color', 'red', 'LineStyle', '--');

% Plot the median as a VTV blue line
h7 = plot(time_period, debt_path_success(:, 5), 'Color',' #ff6875', 'LineWidth', 2,'marker','*'); % 50th percentile, black line
% Plot the median as a black line
%plot(time_period, debt_path_success(:, 5), 'k', 'LineWidth', 2,'marker','*'); % 50th percentile, black line

% Set y-axis limits
%ylim([55 85]);

% convert number to string with space
formatted_simulations = formatWithSpaces(nbr_sim_paths);

if language == 1 % English
    % Customize plot
    xlabel('Year','FontName', 'Calibri');
    ylabel('Debt-to-GDP ratio, %','FontName', 'Calibri');

    if solution_columnZ ~= 0
        titleStr = sprintf('%d%% of the Debt Paths Declining – Annual Adjustment %.2f pp.', plausibility*10, optimal_stoch_a);
        title(titleStr,'FontName', 'Calibri');
    else
        title(sprintf('%d%% of the Debt Paths Declining – Adjustment Plan Not Found', plausibility*10));
    end


    subtitle(sprintf('%d-year adjustment period (%s simulations)', adj_periods, formatted_simulations));
    legend([h7, h5, h6, h4, h3, h2, h1],{'Median','Adjustment plan begins',...
        'Adjustment plan ends','Decile 4–6','Decile 3–7 ',...
        'Decile 2–8','Decile 1–9'}, ...
        'Location','eastoutside','FontName', 'Calibri');
    legend('boxoff')
    %grid on;

    % Set the figure size
    set(gcf, 'Position', [200 200 1.3*FigSize FigSize]);
    
if plausibility == 7
    print(sprintf('DebtFanChart70_%s', methodName), '-dpng', '-r300', '-cmyk');
elseif plausibility == 8
    print(sprintf('DebtFanChart80_%s', methodName), '-dpng', '-r300', '-cmyk');
elseif plausibility == 9
    print(sprintf('DebtFanChart90_%s', methodName), '-dpng', '-r300', '-cmyk');
end


elseif language == 2 % Suomi
    % Customize plot
    xlabel('Vuosi','FontName', 'Calibri');
    ylabel('Velkasuhde, %','FontName', 'Calibri');

    if solution_columnZ ~= 0
        titleStr = sprintf('%d%% velkaurista laskevia – vuosittainen sopeutus %.2f %%-yksikköä', plausibility*10, optimal_stoch_a);
        title(titleStr,'FontName', 'Calibri');
    else
        title(sprintf('%d%% velkaurista laskevia – vuosittaista sopeutusta ei löydy', plausibility*10));
    end
    
    subtitle(sprintf('%d vuoden sopeutusjakso (%s simulaatiota)', adj_periods, formatted_simulations));
    legend([h7, h5, h6, h4, h3, h2, h1],{'Mediaani','Sopeutusjakso alkaa',...
        'Sopeutusjakso päättyy','4. - 6. desiili','3. - 7. desiili',...
        '2. - 8. desiili','1. - 9. desiili'}, ...
        'Location','eastoutside','FontName', 'Calibri');
    legend('boxoff')
    %grid on;

    % Set the figure size
    set(gcf, 'Position', [200 200 1.3*FigSize FigSize]);
    
if plausibility == 7
    print(sprintf('Velkaviuhka70_%s', methodName), '-dpng', '-r300', '-cmyk');
elseif plausibility == 8
    print(sprintf('Velkaviuhka80_%s', methodName), '-dpng', '-r300', '-cmyk');
elseif plausibility == 9
    print(sprintf('Velkaviuhka90_%s', methodName), '-dpng', '-r300', '-cmyk');
end

elseif language == 3 % Swedish
    % Customize plot
    xlabel('År','FontName', 'Calibri');
    ylabel('Skuldkvot, % ','FontName', 'Calibri');

    if solution_columnZ ~= 0
        titleStr = sprintf('%d%% av skuldbanorna minskar – årlig anpassning %.2f procentenheter', plausibility*10, optimal_stoch_a);
        title(titleStr,'FontName', 'Calibri');
    else
        title(sprintf('%d%% velkaurista laskevia - vuosittaista sopeutusta ei löydy', plausibility*10));
    end
    
    subtitle(sprintf('%d års anpassningsperiod (%s simuleringar)', adj_periods, formatted_simulations));
    legend([h7, h5, h6, h4, h3, h2, h1],{'Median','Anpassningsperioden börjar',...
        'Anpassningsperioden avslutas','Decil 4–6','Decil 3–7',...
        'Decil 2–8','Decil 1–9'}, ...
        'Location','eastoutside','FontName', 'Calibri');
    legend('boxoff')
    %grid on;

    % Set the figure size
    set(gcf, 'Position', [200 200 1.3*FigSize FigSize]);
    
if plausibility == 7
    print(sprintf('Skuldfigur70_%s', methodName), '-dpng', '-r300', '-cmyk');
elseif plausibility == 8
    print(sprintf('Skuldfigur80_%s', methodName), '-dpng', '-r300', '-cmyk');
elseif plausibility == 9
    print(sprintf('Skuldfigur90_%s', methodName), '-dpng', '-r300', '-cmyk');
    
end
end

if saveFlag == 1
    save('LastStochModel.mat','-v7.3','-nocompression')
end
end
