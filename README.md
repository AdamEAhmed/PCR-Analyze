# PCR-Analyze

This is an R script for use with quantitative PCR, or real-time PCR, (qPCR).

## Instructions:

***All names entered are case-sensistive***


1. Enter raw CT values into the "RAW" tab of the template spreadsheet, well names are mandatory.  Only numeric data can go here.  If there is a "late call" or e.g. 35, these MUST be left BLANK.


2. Enter the appropriate information into the "KEY" tab, well names and dates are mandatory. Any row left BLANK in the RAW sheet MUST be filled with NA in the KEY sheet under the "date",  "gene", "treatment", and "technical replicate" columns.  The well columns can be left as the actual well number.
	- Label technical replicates accordingly, delta-Ct 	values will be averaged by combining technical replicates first, creating one biological replicate (then biological replicates will be averaged together).


3. In the "GENES" tab, gene 1 is the housekeeping gene and gene 2 is the gene being tested, then in the third column, enter which column number in the "GROUPS" tab references the proper control conditions. (NB - column A = 1, column B = 2, column C = 3, etc.)


4. Enter each control group name into the second row in the "GROUPS" tab, then each of the treatment groups (will be sorted in that order) under the respective control groups.
	- If there are any groups that you would like to be ignored before calculating delta-delta CT values, in last column (after last control/treatment pairing column) row 1, enter "Ignored" (CASE SENSITIVE), and list the ignored groups down this column.
	- Unrelated control groups must have unique names as well as the treatment group's names must also be unique.


5. Save the files by clicking File > Save as, choose the name and location to save the file to, then click on "Save as type:  " and select "Text (Tab delimited) (*.txt)" and click save.


6. Open R.exe

7. Click FIle>Open script... and select the PCR script.

8. If you would not like to check for outliers, at line 24 set "omit" to 0.

9. Click Edit> Run all

10. It will ask to select all of the files that were just saved (raw>key>genes>groups), then at the end it will ask for a folder location to save all the plots and data sets that were produced.
