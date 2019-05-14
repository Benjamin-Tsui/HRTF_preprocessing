# HRTF_preprocessing

### Make sure your have Mapping Toolbox installed
This is a Matlab tool that could perform the required HRTF preparation process for machine learning in a fast and robust way. It could check common errors in SOFA files (e.g. missing data, abnormal measurement distribution), modifies SOFA files (e.g. remove abnormal data), finds common measurement angles within a range, finds the best normalise attributes (HRIR length, sampling frequency), normalises and save the modified measurements as new SOFA files.

This tool consists of different functional building blocks. The main function of this toolbox is preprocess_SOFA.m. The input variables consist of: SOFA files or related folders (input_files), range of tolerance in angle matching (angle_range), measurement distance error tolerance (dist_range) and the name of output folder(output_dir). The function will automatically find the best matching angles, normalise with the proper attributes and save them in a new folder as dictated by the user.

### Main function: preprocess_SOFA.m

Input (with examples): 
1. input = {'ITA_HRTF_Database/SOFA'};     % folder that contains sofa file
    or
    input = {'hrtf b_nh15.sofa'};     % sofa file

 2. angle_range = 0.01;     % set angle matching tolerance (in degree)
    (default angle_range = 0)

 3. dist_range = 0.1;     % set distance tolerance (in meters)  
    (default dist_range = 0)

 4. output_dir = 'normalised_SOFA/';  % set output folder location
     (Optional but highly recommanded)

  Optional input:
 5. plot_trig = 1;     % plot matched angles trigger (1 or 0, default = 1 (on))
 
 6. target_length = 512;     % normalise hrir target length

 7. target_fs = 44100;     % normalise hrir target length

 8. target_fs = 0.99;     % normalise hrir maximum amplitude


Output:
 1. normalised SOFA files with matched angles in the designated folder with an output_file_log list out the angles and modification of each SOFA file
 
 2. summary      % summarise key variables in the process
 (a copy will save in the folder in .mat file)
 
* Detail descriptions are in the script

### Pipeline:
1. Check SOFA file.
2. Organise and fix abnormal SOFA files (may require user’s input to decide what to keep if there is an abnormal SOFA file),
3. Find SOFA files that can represent other files with a similar angular distribution (to improve speed).
4. Find angular match in ‘unique’ SOFA files (may require user’s adjustments if match not found).
5. Find suggested normalisation attributes.
6. Extract measurements, normalise and save as new SOFA files to the desired output folder.
7. Output summary of the new files. (new file names, remaining measurement angles and notes about the file if there is any modification)


### The function signal flow:

![image](https://user-images.githubusercontent.com/25059141/42376052-ab8a15c0-8115-11e8-9c94-8dbdaee9e192.png)


### Using the toolbox:

The suggested input parameters in the main function preprocess_SOFA.m includes SOFA files or related folders (input_files), range of tolerance in angle matching (angle_range), measurement distance error tolerance (dist_range) and the name of output folder (output_dir). Before going through the preparation process, the function will first check all the input attributes, especially whether the appointed output folder already exists. If the folder exists, users can decide whether to overwrite it.

During the organising and fixing process, a folder named good_SOFA will be created to group all the good and repaired SOFA files in one place. If any problematic SOFA files need modification, it will request further instructions in the command window. Otherwise, the function will pause after grouping all the SOFA files. 

By pressing any key to continue, the function will start finding matched angles inside the good_SOFA folder. To avoid spending time on processing the files with similar angular distributions, the function will first run a quick check to find out which SOFA files share the same distribution, then pick one file to represent each distribution in the matching process. During the angle matching process, there is a chance that the function may not find any match, which then requires users to set a wider tolerance range or remove one of the files. After the matching process, the function will pause again and plot matched angles of each represented data set on a graph with a uniform distance for better comparison (unless the user has suppressed the plot).

By pressing any key to continue, the function will normalise and save the processed HRTFs as new SOFA files. The function will find the default recommended attribute (maximum HRIR length and minimum sampling frequency), then normalise and save the files at the appointed output folder.






# Plot HRTF measurements in 3D

The plot_3d_angle.m is a function allows users to plot all HRTF measurement point in three dimensional. It is similar to the SOFAplotGeometry in the original SOFA function but with more flexibility. Firstly, the input besides accepting pre-loaded SOFA struct, this function also accept directly input SOFA files without loading it, or even just the azimuth angle, elevation angle and distance. Secondly, this function has the flexibility to plot different measurements on the same plot with different colour or markers. It helps in comparing the differences between different measurements. Finally, the marker's locations on the graph are displayed in polar coordinates instead of the default Cartesian coordinates, by selecting the marker on the graph, it will show the polar coordinates of the selected point in three dimensions.

![image](https://user-images.githubusercontent.com/25059141/42374728-fdf9cfa8-8110-11e8-99a9-1eeebe8ac973.png)


### Function: plot_3d_angles

Input (3 option):

 Option 1:
   azimuth angle, elevation angle, distance
   e.g. plot_3d_angles(-45, 30, 1.5)

 Option 2:
   loaded sofa file in struct, row number (optional)
   if 2nd input is empty, it will plot all angles
   e.g. hrtf = SOFAload('irc_1007.sofa'); 
        then
        plot_3d_angles(hrtf) % plot all angles
        or
        plot_3d_angles(hrtf, 1)      % plot row 1 in the sofa file
        or
        plot_3d_angles(hrtf, [1 35 60])      % plot row 1, 35 and 60 in the sofa file
  
 Option 3 (similar to option 2):
   name of the sofa file in struct, row number (optional)
   if 2nd input is empty, it will plot all angles
   e.g. plot_3d_angles('irc_1007.sofa')      % plot all angles
        or
        plot_3d_angles('irc_1007.sofa', 1)      % plot row 1 in the sofa file
        or
        plot_3d_angles('irc_1007.sofa', [1 35 60])      % plot row 1, 35 and 60 in the sofa file


 Set Properties:
   different scatter properties can be set after the input data
   tested with 'Marker', 'MarkerEdgeColor' and 'Legend'
 
