# UniFrac workshop

This workshop is about how to use UniFrac in R. It was presented at the Reid/Burton/Gloor lab Data Club on February 29, 2016. A slide deck was presented based on the draft of my [paper](expanding_the_unifrac_toolbox.pdf), submitted to the [Great Lakes Bioinformatics and the Canadian Computational Biology Conference 2016](https://www.iscb.org/glbioccbc2016). This README was then gone through to demonstrate the use of the enclosed [UniFrac R script](UniFrac.r).

## UniFrac

UniFrac is the distance between two microbiome samples. To calculate UniFrac, you need a phylogenetic tree, and a table of counts per bacterial taxa (or operational taxonomic units) per sample.

### Original (Unweighted) UniFrac

[Lozupone et al first publicized unweighted UniFrac in 2005](http://aem.asm.org/content/71/12/8228.full).

![Unweighted UniFrac](images_for_README/unifrac.png)

To calculate unweighted UniFrac, you take the branches of the phylogenetic tree that are not shared between the two samples and and divide that by the total branch lengths. ** Note that you should rarefy your data to the same sequencing depth for unweighted UnIFrac **

### Weighted UniFrac

In 2007, Lozupone et al published an [exciting modification on unweighted UniFrac](http://aem.asm.org/content/73/5/1576.long), where the branches of the phylogenetic tree are weighted by the difference in taxa abundance.

Weighted UniFrac is calculated as follows:

![Weighted UniFrac](images_for_README/weighted.png)

### Information UniFrac

Information UniFrac weights tree branches by uncertainty. A 50/50 composition is more uncertain than a 90/10 composition. This incorporates the taxa abundance evenness. We introduce Information UniFrac in our [conference paper](expanding_the_unifrac_toolbox.pdf).

The weighting for Information Unifrac is calculated as follows, where `p` is the proportion of the taxa, and the logarithm is base 2:

```
p*log(p)
```

Here's a depiction of unweighted (black line), weighted (red line), and information (blue line) UniFrac weightings. The x-axis is the taxa abundance, and the y-axis is the corresponding weight.

![Weights](images_for_README/weights.png)

### Ratio UniFrac

Ratio UniFrac is weighted by taxa abundance, except that the taxa abundance are divided by the geometric mean of all the taxa abundances in the sample. The geometric mean serves as a baseline taxa abundance. We introduce Ratio UniFrac in our [conference paper](expanding_the_unifrac_toolbox.pdf).

Ratio UniFrac is calculated as follows:

![Ratio UniFrac](images_for_README/ratio.png)

## When to use different types of UniFrac

You should use all of them!

* Unweighted UniFrac is good at showing you when you have low level trends, usually indicative of some sort of contaminant. See the barcode example in the slides, where the samples separate according to which row of the 96 well plate they were on.
 * Don't forget to rarefy your data to the same sequencing depth for unweighted UniFrac!
* Weighted UniFrac is the classic tool used in a lot of microbiome research. If you use any of the below methods, you should also use weighted for comparison. Weighted UniFrac shows you separation proportional to abundance differences for taxa between samples.
* Information UniFrac takes into account the abundance evenness of the whole sample, and can differentiate some outliers missed by classically weighted UniFrac.
* Ratio UniFrac takes into account the baseline abundance of taxa, and can also differentiate some outliers missed by classically weighted UniFrac.

If you have a clear difference between your groups, all the different UniFrac methods will show it. If you have a very small difference between groups, it is possible that you will get misleading results with Unweighted Unifrac, which can vary with rarefaction instance. If you have a middling difference, or outliers, it's good to examine the results from each type of UniFrac and make conclusions based on the way the tool works.

## Using the scripts

### Downloading

Here are three different ways to put these files on your computer.

#### Downloading from the web

Click the Download ZIP button, unzip the folder, and now you have all the files!

#### Downloading using GitHub on the command line

Alternatively, you can use GitHub desktop. I've never used it, but I've heard it's really great and user friendly.

Here is how you can clone my files onto your computer using your command line/terminal. Windows users will need to do this using [Git Bash](https://git-for-windows.github.io/). Change to the directory you want to work in, and run the following command:

```
git clone https://github.com/ruthgrace/ruth_unifrac_workshop.git
```

### Installing

Running the UniFrac script requires [an installation of R](https://www.r-project.org/). You will also need the packages phangorn, zCompositions, and vegan, which you can install from CRAN as follows inside R. Pick a CRAN mirror when prompted (on my computer sometimes it takes a minute for this dialog to pop up).

```
install.packages("phangorn")
install.packages("zCompositions")
install.packages("vegan")
```

### Running

#### Example code

You can run the example code all at once in R like so (make sure your working directory is the same one as the example script, or else include the full file path rather than just the file name):

```
source("example_script.r")
```

You can also run it line by line by pasting the lines in R, and using the commands in the troubleshooting section to see what's going on with the different data structures.

#### Customizing the example code

I have included one example of customizing `example_script.r` in the form of `tongue_tongue_script.r`, which uses a subset of tongue samples from the [Human Microbiome Project](http://hmpdacc.org/). The script `get_comparison_data_from_hmp_set.r` extracts the samples from the full HMP [count table](http://downloads.hmpdacc.org/data/HMQCP/otu_table_psn_v35.txt.gz) and [phylogenetic tree](http://downloads.hmpdacc.org/data/HMQCP/rep_set_v35.tre.gz).

The `tongue_tongue_script.r` includes the following optimizations/changes.

* Filters out all operational taxonomic units (OTUs) which are not at least 1% abundant in at least 1 sample (otherwise there are thousands and thousands of OTUs)
 * Also removes the removed OTUs from the phylogenetic tree
* Calculates all the distance matrices (with the exception of Unweighted UniFrac, which is done separately as it requires rarefied data) using this command:
```
all_distance_matrices <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)
```
* Removes all the taxonomy groups code. That code is useful if you want to highlight samples which are more than 50% one taxa (as is the case in the default vaginal data set, where several samples are dominated by different _Lactobacillus_)
* Changes the names of the output files to start with `tongue_tongue_` to differentiate them from the example output. You can also direct the output to a different folder.

### Troubleshooting

When you're getting an error

* Look at your data
 * Sometimes a previous step messed up and now your data is all NA or NaN
 * Sometimes you've run a filter that's too aggressive and now most of your data is gone
* Run things line by line, checking the results as you go along
 * str(myDataFrame) will show you what's inside your R things
 * head(myDataFrame) and tail(myDataFrame) can be used to check the beginning and the end
 * names(myDataFrame) will show you what the column and row names are, so you can make sure it makes sense
 * summary(myVector) will show you where the minimum, maximum, median, mean, and first and 3rd quartile of your data are.

When it's taking forever

* Make sure that all the taxa that you have in your count table and your phylogenetic tree are necessary for your analysis.
 * The scripts work well if you have a couple hundred taxa, but if you have a couple thousand, you may want to filter out the taxa that are less than 1% abundant in every sample.

If you're really super stuck, make an issue on GitHub (Issues tab > New Issue), and I'll take a look when I have time. Good luck!

Ruth