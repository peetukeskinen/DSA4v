% Function to project the debt for a given a
function [debt_out,g_out,drgdp_out,iir_out,pb_out,spb_out] = project_debt4v(scenario, a, iir, potgdp, og, ...
    epsilon,m,dcoa,sfa,inflation,debt_initial, spb,rgdp,...
    debt_st,debt_total,debt_ltn,debt_lt,i_st,i_lt,share_lt_maturing,...
    PB,SF,pb,stoch_sim,g_sim,pb_sim,iir_sim)

        % THIS FUNCTION CALCULATES PROJECTED DEBT PATH GIVEN INPUTS %
        
        % SELECT SCENARIO TO RUN
        % 1) ADJUSTMENT
        % 2) LOWER SPB
        % 3) ADVERSE r-g
        % 4) FINANCIAL STRESS
        
        % STOCHASTIC SIMULATIONS
        % If no simulations required, set stoch_sim = 1
        % and set g_sim,pb_sim,iir_sim to zero
        % 1) NORMAL SHOCKS
        % 2) BLOCK BOOTSTRAP
        % 
        
        % FUNCTION OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % debt: debt-to-gdp projection
        % g: nominal growth of gdp
        % drgdp: real gdp growth
        % iir_out: implicit interest rate output
        % pb_out: primary balance output
        % spb_out: structural primary balance output
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % FUNCTION INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scenario: select scenario
        % stoch_sim: select stochastic simulation
        % a: adjustment per period (in spb terms)
        % iir: implicit interest rate (percent)
        % potgdp: potential output (in MRD euros)
        % og: output gap (percent of potGDP)
        % epsilon: elasticity of budget balance
        % phi: fiscal multiplier on impact
        % dcoa: change in cost of ageing relative to end of adjustment
        % sfa: stock-flow adjustment
        % inflation: price inflation
        % debt_initial: initial level of debt-to-gdp in t
        % spb: COM forecasted spb up to period t+2 (percent of nom GDP)
        % rgdp: COM forecasted real gdp level up to period t+2 (in MRD euros)
        % debt_st: short term debt (in MRD euros)
        % debt_total: total debt (in MRD euros)
        % debt_ltn: new long term debt (in MRD euros)
        % debt_lt: long term debt (in MRD euros)
        % i_st: short term market interest rate
        % i_lt: long term market interset rate
        % share_lt_maturing: share of long term debt maturing yearly
        % PB: primary balance (in MRD euros)
        % SF: stock flow adjustment (in MRD euros)
        % pb: COM forecasted pb up to period t+2 (percent of nom GDP)
        % g_sim: simulated nominal gdp growth shock/draw for stochastic analysis
        % pb_sim: simulated primary balance shock/draw for stochastic analysis
        % iir_sim: simulated implicit interest rate shock/draw for stochastic analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%% Housekeeping        
       
    % Create shells
    periods = length(og); % total number of periods between 2022-2038 (2041)
    drgdp = NaN(periods, 1); % real gdp growth
    g = NaN(periods, 1);     % nominal gdp growth
    debt = NaN(periods, 1);  % debt-to-gdp ratio
    alpha = NaN(periods, 1); % share of short term debt / total debt
    beta = NaN(periods, 1);  % share of new long term debt / total lt debt
    iir_lt = NaN(periods, 1);% long term implied interest rate
    interest_st = NaN(periods, 1);     % interest payment short
    amortization_st = NaN(periods, 1);     
    interest_lt = NaN(periods, 1);     % interest payment long  
    amortization_lt = NaN(periods, 1);
    interest = NaN(periods, 1);        % total interest payments
    amortization = NaN(periods, 1);
    gfn = NaN(periods, 1);             % gross financing needs 
    spb_out = spb;           % shell for spb output vector
   
    % Parameters
    SPB_shock = 0.25;        % lower SPB shock in the first step
    SPB_shock2 = 2*0.25;     % lower SPB shock in the second step
    interest_rate_shock = 0.5; % interest rate shock used in scenario 3&4
    debt(1) = debt_initial;  % initial debt-to-gdp ratio in 2022
    share_st = debt_st(1)/debt_total(1); % share of short-term debt in 2022
    
    ameco_end_t = 3; % last period of ameco forecast
    plan_start_t = ameco_end_t + 1; % adjustment plan start period
    plan_end_t = ameco_end_t + 4; % adjustment plan end period
    monitor_start_t = plan_end_t + 1; % monitoring period start
    
%% Main projection loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = 2:periods % Do projections from period 2 on
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate values for 2023-2024 (pre-adjustment period) %%%%%%%%%
        if t <= ameco_end_t
            rgdp(t) = (1 + og(t)/100) * potgdp(t); %real gdp level
            pb(t) = spb(t) + epsilon*(og(t)); %primary balance
            spb_out(t) = spb(t); %save current period t spb
        
        % calculate values for 2025-2028 (2031) (adjustment period)
        elseif t >= plan_start_t  && t <= plan_end_t 
            pb(t) = (t-ameco_end_t)*a + spb(ameco_end_t) + epsilon*(og(t)-m(t)*a); %primary balance
            spb_out(t) = (t-ameco_end_t)*a + spb(ameco_end_t); %save current period t spb 
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        % Calculate values for 2029-2038 (2041) (post-adjustment period) %
        % depends on the selected scenario
        
            % ADJUSTMENT scenario
        elseif t >= monitor_start_t &&  scenario == 1 % PB post-adjustment period
            % age-adjusted primary balance
            pb(t) = (plan_start_t)*a + spb(ameco_end_t) + epsilon*(og(t)-m(t)*a) + dcoa(t);
            spb_out(t) = (plan_start_t)*a + spb(ameco_end_t); %save current period t spb

            
            % LOWER SPB scenario (step 1)
        elseif t == monitor_start_t && scenario == 2 % 0.25% lower SPB*
            % age-adjusted primary balance
            pb(t) = (plan_start_t)*a + spb(ameco_end_t) -SPB_shock + epsilon*(og(t)-m(t)*a) + dcoa(t);
            spb_out(t) = (plan_start_t)*a + spb(ameco_end_t) -SPB_shock; %save current period t spb
            
            % LOWER SPB scenario (step 2)
        elseif t > monitor_start_t  && scenario == 2 % permanently 0.5% lower SPB*
            %primary balance
            pb(t) = (plan_start_t)*a + spb(ameco_end_t)-SPB_shock2 + epsilon*(og(t)-m(t)*a) + dcoa(t);
            spb_out(t) =  (plan_start_t)*a + spb(ameco_end_t) -SPB_shock2; %save current period t spb
            
            
            % ADVERSE r-g scenario      
        elseif t >= monitor_start_t && scenario == 3 % permanently +0.5% i_st/i_lt
            i_lt(t) = i_lt(t) + interest_rate_shock; % higher short term interest rate 
            i_st(t) = i_st(t) + interest_rate_shock; % higher long term interest rate
            
            % age-adjusted primary balance
            pb(t) = (plan_start_t)*a + spb(ameco_end_t) + epsilon*(og(t)-m(t)*a) + dcoa(t);
            spb_out(t) = (plan_start_t)*a + spb(ameco_end_t); %save current period t spb
            
            
            % FINANCIAL STRESS scenario (step 1)  % temporarely +1% i_st/i_lt
        elseif t == monitor_start_t && scenario == 4 
            i_lt(t) = i_lt(t) + 2*interest_rate_shock; % higher short term interest rate 
            i_st(t) = i_st(t) + 2*interest_rate_shock; % higher long term interest rate
            % age-adjusted primary balance
            pb(t) = (plan_start_t)*a + spb(ameco_end_t) + epsilon*(og(t)-m(t)*a) + dcoa(t);
            spb_out(t) = (plan_start_t)*a + spb(ameco_end_t); %save current period t spb
            
            % FINANCIAL STRESS scenario (step 2)  
        elseif t > monitor_start_t && scenario == 4
            
            % age-adjusted primary balance
            pb(t) = (plan_start_t)*a + spb(ameco_end_t) + epsilon*(og(t)-m(t)*a) + dcoa(t);
            spb_out(t) = (plan_start_t)*a + spb(ameco_end_t); %save current period t spb
        
        else
            disp('No scenario selected.');
            return;
        end
        
        %% Calculate the shares of short-term and long-term debt in total debt
        alpha(t-1) = debt_st(t-1) / debt_total(t-1);
        beta(t-1) = debt_ltn(t-1) / debt_lt(t-1);

        % Use implied iir_lt in the short term and back out iir_lt (year<2025)
        if t <= ameco_end_t
            iir_lt(t) = (iir(t) - alpha(t-1) * i_st(t)) / (1 - alpha(t-1));
        else
            % Use DSM 2023 Annex A3 formulation after (year >= 2025)
            iir_lt(t) = beta(t-1) * i_lt(t) + (1 - beta(t-1)) * iir_lt(t-1);
            iir(t) = alpha(t-1) * i_st(t) + (1 - alpha(t-1)) * iir_lt(t);
        end
            
        %% Interest and amortization payments
        % Short-term debt
        interest_st(t) = debt_st(t-1) * i_st(t-1) / 100; % interest payments on newly issued short-term debt
        amortization_st(t) = debt_st(t-1); % amortization payments on short-term debt share in last year's gross financing needs

        % Long-term debt
        interest_lt(t) = iir_lt(t) / 100 * debt_lt(t-1); % lt interest is t-1 lt debt times implicit lt interest rate
        amortization_lt(t) = share_lt_maturing(t) * debt_lt(t-1); % lt amortization based on maturing share and inst debt

        % Aggregate interest and amortization payments
        interest(t) = interest_st(t) + interest_lt(t); % interest payments on newly issued debt and outstanding legacy debt
        amortization(t) = amortization_st(t) + amortization_lt(t); % amortization of newly issued st and lt debt

        % Gross financing needs based on old and new financing costs and primary balance
        gfn(t) = max([interest(t) + amortization(t) - PB(t) + SF(t), 0]);

        % New debt stock based on new issuance and debt stock last year
        % Total debt
        debt_total(t) = max([debt_total(t-1) - amortization(t) + gfn(t), 0]);

        % Distribution of short-term and long-term debt in financing needs
        D_stn_theoretical = share_st * debt_total(t); % st debt to keep share equal to share_st
        D_ltn_theoretical = (1 - share_st) * debt_total(t) - debt_lt(t-1) + amortization_lt(t); % lt debt to keep share equal to 1 - share_st
        share_st_issuance = D_stn_theoretical / (D_stn_theoretical + D_ltn_theoretical); % share of st in gfn

        % Calculate short-term and long-term debt issuance
        debt_st(t) = share_st_issuance * gfn(t);
        debt_ltn(t) = (1 - share_st_issuance) * gfn(t);
        debt_lt(t) = debt_lt(t-1) - amortization_lt(t) + debt_ltn(t);
        
        %% Calculate implied real gdp level based on og, m, a, potgdp
        % og and a are divided by 100 to get the decimal form
        % real gdp level is affected by adjustment a times the multiplier m
        rgdp(t) = (1 + og(t)/100 - m(t)*a/100) * potgdp(t); % real gdp level 
        
        % Get nominal gdp growth = real gdp growth + inflation
        drgdp(t) =  ((rgdp(t) - rgdp(t-1)) / rgdp(t-1)); % real gdp growth
        
        % adjustment or lower SPB
        g(t) = ((rgdp(t) - rgdp(t-1)) / rgdp(t-1)) + inflation(t); % nominal gdp growth
        
        if t >= monitor_start_t && scenario == 3 % adverse r-g
            g(t) = g(t) - 0.005; % 0.5 %-points lower nominal gdp growth
        end
        
        %% Project debt using debt dynamic equation  
        
        if stoch_sim == 1 || stoch_sim == 2
        % define debt drivers with normal shocks
        g(t) = g(t) + g_sim(t);
        pb(t) = pb(t) + pb_sim(t);
        iir(t) = iir(t) + iir_sim(t);
        
        % BOOTSTRAP UNCOMMENTED
%         elseif stoch_sim == 2 && t >= 8 && t <= 12
%         % define debt drivers with block-bootstrap shocks
%         g(t) = g_sim(t);
%         pb(t) = pb_sim(t);
%         iir(t) = iir_sim(t);
%         
%         elseif stoch_sim == 2 && (t < 8 || t > 12)
%         % outside stochastic simulation period
%         g(t) = g(t);
%         pb(t) = pb(t);
%         iir(t) = iir(t);
        
        else
            disp('Stochastic simulation method not set.');
            return;
        end
        
        % equation 3 in the report, p. 57)
        debt(t) = debt(t-1) * ((1 + iir(t)/100) / (1 + g(t))) - pb(t) + sfa(t);
        
 
    end
    
    %% Save output variables
    debt_out = debt;
    g_out = g;
    drgdp_out = drgdp;
    iir_out = iir;
    pb_out = pb;
    
    
end