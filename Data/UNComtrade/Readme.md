This folder contains the data needed to construct a global trade network.  

The primary data set is "UNComtrade.csv", which contains information on 
commodities flows (i.e., trade in goods, not services) involving 179 reporter 
countries and 245 partner countries for 2016.
Publicly available data from https://comtradeplus.un.org/ (free, requires login)
were merged with region information to produce the included dataset. 

Included variables are: 

- Year
- ReporterISO
- ReporterDesc
- ReporterRegion
- ReporterSubRegion
- PartnerISO
- PartnerDesc
- PartnerRegion
- PartnerSubregion
- FlowDesc: import or export
- Primary value: value of the trade flow

Additional information on methods available here: 
https://comtrade.un.org/data/MethodologyGuideforComtradePlus.pdf and
https://unstats.un.org/wiki/display/comtrade/New+Comtrade+FAQ+for+First+Time+Users


