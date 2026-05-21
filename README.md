# Dual-Dunnett-Barplot
A MATLAB program to handle multiple dunnetts test across dimensions of a square matrix of barplot data (genotypes * conditions)

This readme was last edited: 2026-05-21


<img width="480" height="480" alt="ExampleData" src="https://github.com/user-attachments/assets/ae785b6f-5c7c-406c-8dfc-13aefb5d547a" />


# Introduction: The Why
Reading academic journal articles, I enjoy seeing clean figures that show the data in one panel, with annotations for statistical significance below the data display. I wanted to create something that captured this format (as it can be produced in GraphPad prism easily, or manually, but less so programmatically in MATLAB. This function is currently usable for a 3x3 set of genotypes and conditions, and provides a compact visualization of the data and statistical comparisons for a dunnett's test multiple comparisons (hence the _current_ name).

# Usage:

A selected snipped from the top of the code shows the main input variables:
```
function dual_dunnett_barplot(data,genotypes,concentrations,cmap,upperYLabel,defaultsettings)
% ...
if nargin<6
    defaultsettings = [
    {[3,1]}, % tiledlayout size 
    {[1,1]}, % upper bar plot tile size
    {[2,1]}, % lower annotations tile size
    {abyss(num_conc)}, % annotations colormap
    ];
end
```
* Inputs:
  * data: 3D matrix of datapoints (dim 1), conditions (dim 2), and genotypes (dim 3)
  * genotypes: cell array of strings corresponding to the size of dim 3 of the data matrix
  * concentrations: cell array of strings corresponding to the size of dim 2 of the data matrix
  * cmap: colormap for bars in the upper plot of the resulting figure, a ```brewermap(...)``` provides a good set of colors for most use
  * upperYLabel: Y-label for the bar plot to have the text auto-populate on the figure
  * defaultsettings: cell array of small matrices defining plot and annotation specific configuration details


This script uses input validation to ensure the user is giving data that will be able to plot sucessfully:
```
% Input Validation
if (length(genotypes)~=size(data,3)) error("Mismatch between size of genotypes array and dim 3 of input data..."); end
if (length(concentrations)~=size(data,2)) error("Mismatch between size of concentrations array and dim 2 of input data..."); end
if (size(cmap,1)~=size(data,2)) error("Mismatch between colormap and num of concentrations in input data..."); end
```

This script requires the following helper function to convert p-values to stars (provided):
```
function stars = p_to_stars(p)
if p <= 0.001,     stars = '***';
elseif p <= 0.01,  stars = '**';
elseif p <= 0.05,  stars = '*';
else,              stars = 'ns';
end
end
```
Also needed is a function to return the Dunnett Multiple Comparisons after 1-way ANOVA (provided):
```
function mc = returnDunnettMultCompareStats(data)
[p,anovatab,stats] = anova1(data(:,:,1),[],"off");
mc = multcompare(stats,"ControlGroup",1,"Approximate",1,"CriticalValueType","dunnett","Display","off","Alpha",0.05/12);
end
```

# Additional info:

* For ```brewermap(...)```: https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
* It is helpful to define the figure and size appropriately first before running the actual plot function. As exampled:
```
figure(1);clf;set(gcf,"Units","inches","Position",[1,1,5,5]);
```


## Roadmap:
1. Add in ability to use a user defined size of input data matrix
2. Fix handling of "ns" and stars so that the text is 'nicely' aligned to the corresponding spreader bars
3. Incorporate options to run different secondary statistical analyses within conditions, but across genotypes >> in case other options are more desirable or are more statistically valid
