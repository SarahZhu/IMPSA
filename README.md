This the implementation code of my master project at University of California, Riverside under my advisor Dr. Eamonn Keogh:  *Implementation and Optimizations of the Matrix Profile Shapelet Algorithm on Time Series Classification*, which is based on lots of previous works  which I have commented in the code. A more general introduction of works on Matrix Profile is here: *https://www.cs.ucr.edu/~eamonn/MatrixProfile.html* .

1. Some datasets we used in our experiment
	* FuzzySine128
	![FuzzySine128](https://github.com/SarahZhu/IMPSA/tree/master/images/FuzzySine128.jpg)
	* StarLightCurves
	![StarLightCurves](https://github.com/SarahZhu/IMPSA/tree/master/images/StarLightCurves.jpg)
2. Comparisons among optimizations and the baseline: total time, test accuracy and time bar graph showed in two major parts.
	* Total time
	![total time](https://github.com/SarahZhu/IMPSA/tree/master/images/Comparison_of_total_time_among_optimization_methods.jpg)
	* Accuracy
	![accuracy](https://github.com/SarahZhu/IMPSA/tree/master/images/Comparison_of_accuracy_among_optimization_methods.jpg) 
	* Total time showed by two major components in the program
	![time bar](https://github.com/SarahZhu/IMPSA/tree/master/images/Comparison_of_total_time_showed_in_2_parts_among_optimization_methods.jpg)

3. Introduction of critical modules:

Module Name | Description | My Contribution
----------- | ----------- | ---------------
main | Matrix Profile Shapelet Algorithm (MPSA) | Fixed the instability issue of occurring NaNs and Infs occasionally in the original code; Fixed small bugs in plotting shapelet.
main2 | Iterative Matrix Profile Shapelet Algorithm (IMPSA) | Plugged MPSA into loops, looping thru all shapelet candidates of all possible lengths, with early abandon strategy added
main3 | IMPSA with output of both train and test accuracy, early abandon strategy emitted to search for candidates of all possible lengths | Changed a few lines to plot train and test accuracy for all shapelet of searched lengths
PreSCRIMP | An existing version of PreSCRIMP on Matrix Profile self-join P<sub>AA</sub> | Added handlers for NaNs and Infs in Matrix Profile for MPSA
PreSCRIMP_joinAB | PreSCRIMP on Matrix Profile join P<sub>AB</sub> | Implemented this version based on the self-join version
prepare_data_tsv | Prepare input time series data in tsv format | Self-developed
MASS_V2 | The MASS method | Fully reference to Mueen's Similarity Search
generate_sin_ts | Generate FuzzySine dataset | Self-developed



4. How to run the code:
	* First please download the time series data set from UCR Time Series Archive (https://www.cs.ucr.edu/~eamonn/time_series_data_2018/ ) 
	* After downloading those data, upzip and put data (tsv format) in the ***data*** folder. There are already some sample datasets in the folder for you to get used to the program.
	* Modify the sublength in the ***main*** program to try MPSA, or:
	* Modify ***ts_len*** as the time series length corresponding to the dataset you are using, and ***step*** in ***main2*** or ***main3*** to changed the interval of your search on the subsequence length of shapelet candidates. The larger the step, the faster you will get the result, but lower the accuracy. Carefully choose the step size to find a balance between speed and accuracy.
	* Note that this version of IMPSA only handles binary classification using only one shapelet to build the decision tree. However, it's trivial to use multiple shapelets for a multi-class classification.

Please feel free to contact me at sarah.zhu@email.ucr.edu if you have any questions on running the code.
