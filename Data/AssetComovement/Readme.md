This folder contains the data needed to construct an asset price co-movement network for a subset of NASDAQ tickers. 

The primary data set is Stockmarket.csv, which contains information on 182 NASDAQ 
stocks included on the index from 2015-2020. 
Publicly available data from https://www.kaggle.com/datasets/paultimothymooney/stock-market-data/data, 
subset for the time indicated, and processed to produce the following variables: 

- Month
- Year
- Symbol
- Volume: mean monthly volume
- MonthlyReturns: simple (arithmetic) monthly returns computed using closing price 

The secondary data set consists of the Fama French 5 factors indices for the same period, 
publicly available here: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html#Research
Included variables: 

- Month
- Year
- Mkt.RF: value-weight return of all CRSP firms 
- SMB: avg return on the nine small stock portfolios - avg return on the nine big stock portfolios
- HML: avg return on the two value portfolios minus the avg return on the two growth portfolios
- RMW: avg return on the two robust operating profitability portfolios minus the avg return on the two weak operating profitability portfolios
- CMA: avg return on the two conservative investment portfolios minus the avg return on the two aggressive investment portfolios
- RF: Treasury bill rate

Additional information on variables is available here: 
https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/f-f_5_factors_2x3.html

Many analyses are possible using the above data. One possible approach is to:
- Compute idiosyncratic returns for each stock by regressing monthly returns on the FF5 as in the Fama French model
- Compute correlation of idiosyncratic returns for all pairs of stocks
- Define a network by thresholding the correlation coefficient
- Extension: repeat the following at a finer time scale to create a time series of networks









