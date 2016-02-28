# UniFrac workshop

This workshop is about how to use UniFrac in R. It was presented at the Reid/Burton/Gloor lab Data Club on February 29, 2016. A slide deck was presented based on the draft of my [paper](expanding_the_unifrac_toolbox.pdf), submitted to the [Great Lakes Bioinformatics and the Canadian Computational Biology Conference 2016](https://www.iscb.org/glbioccbc2016). This README was then gone through to demonstrate the use of the enclosed [UniFrac R script](UniFrac.r).

## UniFrac

### Original (Unweighted) UniFrac

cite lozupone

include image

### Weighted UniFrac

cite lozupone

### Exponent UniFrac

link to plos manuscript

### Ratio UniFrac

link to plos manuscript

## When to use different types of UniFrac

You should use all of them!

* Unweighted UniFrac is good at showing you when you have low level trends, usually indicative of some sort of contaminant. See the barcode example in the slides, where the samples separate according to which row of the 96 well plate they were on.
* Weighted UniFrac is the classic tool used in a lot of microbiome research. If you use any of the below methods, you should also use weighted for comparison. Weighted UniFrac shows you separation proportional to abundance differences for taxa between samples.
* Information UniFrac takes into account the abundance profile of the whole sample, and can differentiate some outliers missed by classically weighted UniFrac.
* Ratio UniFrac takes into account the baseline abundance of taxa, and can also differentiate some outliers missed by classically weighted UniFrac.

If you have a clear difference between your groups, all the different UniFrac methods will show it. If you have a very small difference between groups, it is possible that you will get misleading results with Unweighted Unifrac, which can vary with rarefaction instance. If you have a middling difference, or outliers, it's good to examine the results from each type of UniFrac and make conclusions based on the way the tool works.

## Using the scripts

### Downloading

Here are three different ways to put these files on your computer.

#### Downloading from the web

#### Downloading using GitHub Desktop

#### Downloading using GitHub on the command line

### Installing

Running the UniFrac script requires an installation of R. You will also need the packages phangorn, zCompositions, and vegan, which you can install from CRAN as follows inside R. Pick a CRAN mirror when prompted (on my computer sometimes it takes a minute for this dialog to pop up).

```
install.packages("phangorn")
install.packages("zCompositions")
install.packages("vegan")
```

### Running

#### Example code

#### Customizing the example code

### Troubleshooting

When you're getting an error

* look at your data
 * sometimes a previous step messed up and now your data is all NA or NaN
 * sometimes you've run a filter that's too aggressive and now most of your data is gone
* run things line by line, checking the results as you go along
 * str(myDataFrame) will show you what's inside your R things
 * head(myDataFrame) and tail(myDataFrame) can be used to check the beginning and the end
 * names(myDataFrame) will show you what the column and row names are, so you can make sure it makes sense
 * summary(myVector) will show you where the minimum, maximum, median, mean, and first and 3rd quartile of your data are.
 


