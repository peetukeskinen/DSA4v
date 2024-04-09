function [YearlyShocks] = sumq2y(QuarterlyShocks,m_res_lt,Years)
%   SUM quarterly shocks to yearly shocks given the persistence
%   QuarterlyShocks are shocks generated
%   m_res_lt is the average residual maturity
%   N is the number of 5-year forecasts
%   Years is the number of years in each forecast

weight = 1/m_res_lt; 
Q = 4; % number of quarters in a year
TotQuarters = Years*Q; % total quarters in simulation
TotShocks = length(QuarterlyShocks); % number of generated shocks
N = TotShocks/TotQuarters; % number of 5-year forecasts
YearlyShocks = zeros(N,Years); % shell for aggregated yearly shocks

% reshape 
QuarterlyShocks = reshape(QuarterlyShocks,[TotQuarters N]);

for sim = 1:N % simulations
for year = 1:Years % years

YearlyShocks(sim,year) = year*weight*sum(QuarterlyShocks(1:year*Q,sim));
%%%%%%%%%%%%CHECK THIS%%%%%%%%%%%%%%%%

end
end

