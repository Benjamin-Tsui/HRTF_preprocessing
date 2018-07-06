# HRTF_preprocessing

This is a tool that could perform the required HRTF preparation process for machine learning in a fast and robust way. It could check common errors in SOFA files (e.g. missing data, abnormal measurement distribution), modifies SOFA files (e.g. remove abnormal data), finds common measurement angles within a range, finds the best normalise attributes (HRIR length, sampling frequency), normalises and save the modified measurements as new SOFA files.

The function signal flow:

![image](https://user-images.githubusercontent.com/25059141/42374794-37cadb82-8111-11e8-87bd-37e22e728f38.png)


# Plot HRTF measurements in 3D

This function allows users to plot all HRTF measurement point in three dimensional. It is similar to the SOFAplotGeometry in the original SOFA function but with more flexibility. Firstly, the input besides accepting pre-loaded SOFA struct, this function also accept directly input SOFA files without loading it, or even just the azimuth angle, elevation angle and distance. Secondly, this function has the flexibility to plot different measurements on the same plot with different colour or markers. It helps in comparing the differences between different measurements. Finally, the marker's locations on the graph are displayed in polar coordinates instead of the default Cartesian coordinates, by selecting the marker on the graph, it will show the polar coordinates of the selected point in three dimensions.

![image](https://user-images.githubusercontent.com/25059141/42374728-fdf9cfa8-8110-11e8-99a9-1eeebe8ac973.png)
