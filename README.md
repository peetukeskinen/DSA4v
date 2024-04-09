# Debt Sustainability Analysis (DSA) tool	
### Version 4

This repository contains Debt Sustainability Analysis (DSA) tool used in the debt ratio projections published in the [blog post](https://www.vtv.fi/en/blog/the-length-of-the-adjustment-plan-in-the-reformed-eu-debt-rules-is-of-great-importance-to-finland/). This version 4 updates version 2 published in chapter 3 of [Fiscal Policy Monitoring Report 2023](https://www.vtv.fi/en/publications/fiscal-policy-monitoring-report-2023/) published by NAOF's Fiscal Policy Monitoring Unit. 

The code was successfully run with Windows 10 (64-bit) and MATLAB R2020b. The files needed to run the MATLAB code are:

#### 1)	main function runDsaStochModel4_blog.m,
#### 2a)	helper function project_debt4v.m to project debt paths given yearly adjustment,
#### 2b)	helper function sumq2y.m to sum quarterly shocks to yearly shocks,
#### 2c) helper function formatWithSpaces.m to format numbers in figures,
#### 3)	data file AmecoFinlandData4v_1.xlsx.

The main function runDsaStochModel4_blog.m calculates minimum yearly adjustment for four-year adjustment plan for Finland. Version 4 do not consider (yet) all of the criteria. Based on [Commission Debt Sustainability Monitor 2023](https://economy-finance.ec.europa.eu/document/download/e3a23fba-1402-4cc9-b571-7473b5e7842a_en?filename=ip271_en.pdf), the version includes:

#### DSA-based criteria (deterministic and stochastic scenarios)
#### Debt Sustainability Safeguard ([article 6a](https://www.consilium.europa.eu/media/70386/st06645-re01-en24.pdf))

Resulting output and graphs can be saved in the current folder as .txt and .png files.
 
The main function produces debt projections following the implementation of the European Commission's Debt Sustainability Monitor 2023 published in March 2024. In the current MATLAB implementation, I have extensively followed [analysis](https://www.bruegel.org/working-paper/quantitative-evaluation-european-commissions-fiscal-governance-proposal) and [Python code](https://github.com/lennardwelslau/eu-debt-sustainability-analysis) made by Darvas et al. (2023).

Currently, all four scenarios are available: baseline adjustment scenario, lower structural primary balance (SPB) scenario, Adverse r-g scenario and Financial stress sccenario. Also, stock-flow-adjustment (SFA) method can be chosen to follow Commission assumption, or an alternative assumption based on linearly declining SFA. 

Selection for scenario and SFA method can be made when passing function arguments (see below). Parameter section sets the parameter values to match those values chosen by European Commission. User can modify parameter values to do sensitivity testing.

DSA projections can be made with or without considering debt sustainability safeguard. Plotting can be selected on or off. Number of samples and method used in the stochastic simulations can be set. Also language and saving option is available.

AmecoFinlandData4v_1.xlsx contains all the necessary data to run the main code. Tab COM in the Excel file contains data and tab DataSources describes the data sources and some clarifying definitions. COM data is from Commission 2023 autumn forecast round. STOCH tab includes data used in the stochastic simulations. All data is from publicly available sources expect financial market data which cannot be shared and must be separately downloaded via a Bloomberg Terminal. AMECO data can be accessed here. Variables can be found using AMECO variable codes listed in the Excel file under the tab DataSources (e.g. FIN.1.0.0.0.AYIGD). COM projections for some variables can be found here.

For comments and suggestions, please contact peetu.keskinen@vtv.fi.

## Example

Example to run the main function:
### runDsaStochModel4_blog(0,1,1,1,3,7,1,1,0)

Above command runs the function with COM ZERO ASSUMPTION,ADJUSTMENT,
DEBT SAFEGUARD, PLOTTING, 1000 simulations, 70% plausibility value,
English output, Normal (COM) method and no saving.

The arguments define options when running the function.
These are described below:

    SELECT SFA METHOD:
         0) COM ZERO ASSUMPTION
        -1) ALTERNATIVE ASSUMPTION      
    SELECT SCENARIO:
        1) ADJUSTMENT
        2) LOWER SPB
        3) ADVERSE r-g
        4) FINANCIAL STRESS
    APPLY DEBT SAFEGUARD
        yes=1, no=0
    PLOTTING
        yes=1, no=0
    Stochastic samples (power)
        1 = 10 simulated paths
        2 = 100
        3 = 1000
        4 = 10 0000
        5 = 100 0000
        6 = million
    Plausibility value
        7 = 70%
        8 = 80%
        9 = 90%
    Language for Stochastic plots
        1 = English
        2 = Suomi
        3 = Swedish
    Stochastic Method
        1 = Normal
        2 = Bootstrap
    Save results as -mat (OPTIONAL)
        1 = save


