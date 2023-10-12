# ALAn
Automated Layer Analysis (ALAn) is an image analysis tool that provides a quantiative, unbiased characterization of the architecture of epithelial layers. The tool has been optimized for the analysis of Mardin-Darby Canine Kidney cells grown in culture, fixed, and stained with fluorescent dyes marking F-Actin and DNA. ALAn is written in Python using Jupyter Notebooks for user ease.

ALAn is introduced in Dawney NS*, Cammarota C*, Jia Q, Shipley A, Glichowski JA, Vasandani M, Finegan TM, Bergstralh DT. A novel tool for the unbiased characterization of epithelial monolayer development in culture. Mol Biol Cell. 2023 Apr 1;34(4):ar25  doi: https://doi.org/10.1091/mbc.E22-04-0121 (* denotes equal contribution). The tool was further developed in Cammarota C, Dawney NS, Bellomio PM, Jüng M, Fletcher AG, Finegan TM, Bergstralh DT. The Mechanical Influence of Densification on Initial Epithelial Architecture. bioRxiv 2023.05.07.539758; doi: https://doi.org/10.1101/2023.05.07.539758.

A detailed protocol for the use of ALAn is published: Cammarota, C., Bergstralh, D.T., Finegan, T.M. (2023) Automated Layer Analysis (ALAn): An Image Analysis Tool for the Unbiased Characterization of Mammalian Epithelial Architecture in Culture. Bioprotocols, https://doi.org/10.21769/p2429

This read-me includes a list of each function, with all inputs and outputs as well as a brief mention of where things will be saved or requirements for locations to read in files if necessary. ALAn is to be used with 4D confocal stacks (xyz, multichannel images saved with a 512 by 512 ) in conjunction with a spreadsheet of segmented nuclear positions and sizes (xyz positions, volume in µm3­, and reference positions). Example datasets can be downloaded from: [https://mailmissouri-my.sharepoint.com/:f:/g/personal/tmfhdb_umsystem_edu/EsETAvroQ7BLhZEbFpNE-fEBxR8AdpvZuVMseUwO8Y1KzA?e=lB2y7b](https://mailmissouri-my.sharepoint.com/personal/tmfhdb_umsystem_edu/_layouts/15/onedrive.aspx?web=1&id=%2Fpersonal%2Ftmfhdb%5Fumsystem%5Fedu%2FDocuments%2FExample%20Data%20Sets%20for%20ALAn&FolderCTID=0x012000E6A96FE55D1B60449D2C52A02D571508&view=0)
 
The code block in Cell 3 of the ALAn notebook provides a quick way to read your files into the workspace. To use this portion, you will need to name associated image and spreadsheet files the same, up to the extension (.tif or .csv). By changing the path in line 3 to file path for your images (.tif format) and spreadsheets (.csv format), this block will output three lists (list of names, list of dfs, and list of unshuffled images) which correspond to file names, spreadsheets as pandas dataframes, and images as array objects. The positions of each name/spreadsheet/image in these lists will be identical for your images. Any image without an associated spreadsheet, or vice versa, will be left out of the described lists and can instead be found in the unmatched images or unmatched spreadsheets lists.
 

#### image_shuffle (image)
This function is meant to reorder images such that the array object follows the same order (z, channel, x, y). Images saved as 512 by 512 files can have different orderings depending on how they were saved and this function will standardize all images before moving forward.
INPUTS:
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
OUTPUTS:
o shuffled image: 4D numpy array – Transposed input to force the array into mxnx512x512 shape.
 

#### get_layer_position (df, image)
This function extracts the absolute bounds of the image (for images segmented in imaris this will likely be in terms of the microscope stage used to acquire the image). These bounds should be stored in the fifth column of the DataFrame associated with the image, and are used to interpret the nuclear centroid positions of the segmentation. This function will also return an mx1 numpy array where m is the number of slices in an image. This array is the z-position of each slice in the image shifted such that the bottom-most slize is at z=0.
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
OUTPUTS:
o x_max: float – The maximum x coordinate. 
o y_max: float – The maximum y coordinate. 
o z_max: float – The maximum z coordinate. 
o x_min: float – The minimum x coordinate. 
o y_min: float – The minimum y coordinate. 
o z_min: float – The minimum z coordinate.
o slize_heights: 1D numpy array – Position of each z-slice in an mx1 array. 
 

#### local_density (df, image, box_radius = (x_max+x_min)/2)
This function is meant to return the nuclear density (nuclei x 103/mm2) of a box centered at the center of the image. It defaults to the entire image, but users can optionally input a smaller box to find the nuclear density of a smaller region. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o box_radius: float, optional – A float that should be smaller than half the field of view of the image.   
OUTPUTS:
o number_of_cells_in_box: int – Number of nuclei in a pillar bounded by a square box with radius defined by the user up to the entire image. 
o box_density: float – Density of nuclei in a pillar bounded by a square box with radius defined by the user up to the entire image with units (nuclei x 103/mm2). 
 

#### clear_debris (df, plot = False, save_name = False)
This function removes any segmented object determined not to be a nucleus by comparing the size to the size distribution of nuclei in the image. It has options to plot the objects kept vs. the objects removed as a histogram and the option to save that histogram as a pdf. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o plot: boolean, optional – Default option is to not plot the histogram of segmented objects, but will plot all objects and objects kept if set to True or if save_name is defined.   
o save_name: string, optional – A string specifying the name to save the histogram as if desired. Defaults to False, meaning it will not save.
OUTPUTS:
o df_cleared: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented nuclei from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
 

#### layer_height _actin (df, image, plot = False, save_name = False, actin_channel = 1, section = (False, False, False, False), invert = False)
This function is meant to determine cutoffs for the top and bottom of the layer given the actin distribution. The exact mechanism to determine these cutoffs are described in Dawney and Cammarota et al. 2023. The function outputs the layer top and bottom positions, layer height, position of the peak actin intensity, and a normalized actin intensity array of length mx1. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o plot: boolean, optional – Default option is to not plot the actin distribution, but will if True or if save_name is defined.   
o save_name: string, optional – A string specifying the name to save the actin distribution as if desired. Defaults to False, meaning it will not save.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o section: tuple, optional – A tuple with four elements corresponding to the bottom, top, left, and right positions of a region of interest within the greater field of view. Defaults to False positions which cause the entire field of view to be analyzed.
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o min_layer_height: float – The z position of the slice determined to be the bottom of the cell layer. 
o max_layer_height: float – The z position of the slice determined to be the top of the cell layer. 
o layer_height: float – The height of the cell layer in µm. 
o actin_peak: float – The z position of the most intense actin slice. 
o norm_intensities_actin: 1D numpy array – mx1 array of the normalized actin intensities of each slice. 
 

#### smooth_array (array, number_to_smooth = 3)
This function convolves a uniform window of size number_to_smooth with value 1/number_to_smooth with the given array to output a rolling average. 
INPUTS:
o array: 1D numpy array – An array with shape Sx1.  
o number_to_smooth: odd integer between 3 and S, optional – An odd integer defining the convolution window.  
OUTPUTS:
o new_array: 1D numpy array – An array with shape Sx1, product of a rolling average of the input array. 
 

#### find_shoulders (df, image, plot = False, save_name = False, actin_channel = 1, section = (False, False, False, False), invert = False)
Function to detect whether the actin intensity plot has a ‘shoulder’, the exact mechanism to determine the shoulder is described in Dawney and Cammarota et al. 2023. This function returns the number of peaks and the ratio of the peak heights, or Later-to-Apical shape index.
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o plot: boolean, optional – Default option is to not plot the actin distribution and derivative, but will if True or if save_name is defined.   
o save_name: string, optional – A string specifying the name to save the actin distribution and derivative as if desired. Defaults to False, meaning it will not save.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o section: tuple, optional – A tuple with four elements corresponding to the bottom, top, left, and right positions of a region of interest within the greater field of view. Defaults to False positions which cause the entire field of view to be analyzed.
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o equivalence: float – The Lateral-to-Apical shape index for the input layer.
o number_of_peaks: int – The number of peaks detected by the scipy.signal.find_peaks function used on the derivative the normalized actin intensity plot with a five point rolling average with a prominence of 0.008. 


#### gaussian_fit (x, A, B, C)
This function returns a Gaussian peak as defined by parameters A, B, and C with positions x. 
INPUTS:
o x: 1D numpy array – An array of the input positions to a gaussian curve.   
o A: float – Scalar multiplier to the distribution.   
o B: float – Location of the center of the distribution.  
o C: float – Proxy for the width of the distribution.    
OUTPUTS:
o gaussian_distribution: 1D numpy array – A*exp(-((x-B)/C)2)
 

#### double_gaussian_fit (x, A, B, C, D, E, F)
This function returns a sum of two Gaussian peaks, one defined by parameters A, B, and C and the other defined by parameters D, E, and F, with positions x. 
INPUTS:
o x: 1D numpy array – An array of the input positions to a gaussian curve.   
o A: float – Scalar multiplier to peak 1.   
o B: float – Location of the center of peak 1.   
o C: float – Proxy for the width of peak 1.    
o D: float – Scalar multiplier to peak 2.   
o E: float – Location of the center of peak 2.   
o F: float – Proxy for the width of peak 2.    
OUTPUTS:
o double_gaussian_distribution: 1D numpy array – A*exp(-((x-B)/C)2) + D*exp(-((x-E)/F)2)

 
#### nuclei_distribution (df, image, plot = False, save_name = False, actin_channel = 1, section = (False, False, False, False), invert = False)
This function is meant to determine whether the nuclei are distributed in a single monolayer at the bottom of the cell culture well. Exact determination is described in Dawney and Cammarota et al 2023. This function outputs the position of the peak of the nuclear distribution and the organized vs. disorganized classification. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o plot: boolean, optional – Default option is to not plot the nuclear distribution, but will if True or if save_name is defined.   
o save_name: string, optional – A string specifying the name to save the nuclear distribution as if desired. Defaults to False, meaning it will not save.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o section: tuple, optional – A tuple with four elements corresponding to the bottom, top, left, and right positions of a region of interest within the greater field of view. Defaults to False positions which cause the entire field of view to be analyzed.
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o nuclear_peak: float – Z position of the peak of the nuclear distribution. 
o layer: string – Classifies the layer as ‘Organized’ or ‘Disorganized’
 

#### terrain_map (df, image, color = ‘magma_r’, save_name = False, invert = False)
This function plots all of the nuclei projected onto the XY plane, sized based on nuclear size, and colored based on z position of the nucleus. This allows for visualization of the position of extralayer cells or disorganized regions within the field of view. 
INPUTS:
df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o color: string, optional – Matplotlib color name for the look up table to be used for the nuclear position.
o save_name: string, optional – A string specifying the name to save the terrain map as if desired. Defaults to False, meaning it will not save.
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o plot – The only output is the plot as described above. 
 

#### z_projection_with_cutoffs (df, image, plot = False, save_name = False, actin_channel = 1, section = (False, False, False, False), invert = False)
This function returns the top and bottom actin slices, and can create and XZ plot with the image projected onto the y-axis. This side view is good to test the accuracy of the top and bottom finding function. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o plot: boolean, optional – Default option is to not plot the nuclear distribution, but will if True or if save_name is defined.   
o save_name: string, optional – A string specifying the name to save the nuclear distribution as if desired. Defaults to False, meaning it will not save.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o section: tuple, optional – A tuple with four elements corresponding to the bottom, top, left, and right positions of a region of interest within the greater field of view. Defaults to False positions which cause the entire field of view to be analyzed.
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o min_actin_slice: int – The z slice determined to be the bottom of the cell layer. 
o max_actin_slice: int – The z slice determined to be the top of the cell layer. 
 
#### layer_determination (df, image, actin_channel = 1, section = (False, False, False, False), invert = False)
This function is meant to output the important information for each layer. The first output is the layer determination as described in Dawney and Cammarota et al 2023, and expanded in Cammarota et al. 2023 (preprint). Additional outputs are the number of cells on and in the layer, percentage of extralayer cells, and cell density. For the ‘Disorganized’ classification, the only meaningful quantity is cell density.
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o section: tuple, optional – A tuple with four elements corresponding to the bottom, top, left, and right positions of a region of interest within the greater field of view. Defaults to False positions which cause the entire field of view to be analyzed.
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o layer_classification: string – One of five layer classes, ‘Immature’, ‘Intermediate A’, ‘Intermediate B’, ‘Mature’, and ‘Disorganized’.
o cells_above: int – Integer number of cells above of the cell layer. 
o cells_inside: int – Integer number of cells in the cell layer.
o percentage_above: float – Perceent of cells above the layer.
o cell_density: float – Density of nuclei in a pillar bounded by a square box with radius defined by the user up to the entire image with units (nuclei x 103/mm2).
 
#### sub_classify (image)
This function breaks a full field of view image into the smallest classifiable sections (100px x 100 px), and classifies each of the subsections. In order to make this sub classification, the image removes the final two rows and columns from the 512x512 image to make 25 102x102 images. Each image is run through a layer classification individually to find individual layer heights and other properties. This function plots a grid of 25 squares colored by layer type. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o save_name: string, optional – A string specifying the name to save the sub-classified image as if desired. Defaults to False, meaning it will not save.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o all_sections_analyzed: list – This is a list of 25 five element tuples. Each tuple is the (sub_layer_classification, cells_above_sub_layer, cells_inside_sub_layer, sub_layer_height, super_layer_classification), one for each of the squares in the grid. 
 
#### nuclear_centroid_actin_overlay (df, image, save_name = False, actin_channel = 1)
This function is used to test the segmentation and clear debris functions to ensure that nuclei used by ALAn are reasonable. The nuclear segmentations in the layer are overlayed with the actin channel to view the position of nuclei within each cell. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o save_name: string, optional – A string specifying the name to save the nuclear centroid overlayed onto the actin signal as if desired. Defaults to False, meaning it will not save.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
OUTPUTS:
o plot – The only output is the plot described above.
 
#### xy_segmentation (df, image, plot = False, save_name = False, actin_channel = 1, invert = False)
This function segments the actin of the underlying organized cell layer. This only works for organized layer types. This function sums the slices starting from 3 above the bottom of the layer and goes to three below the top of the layer. Using a series of binarization, opening and closing, and other skimage features, the cells are segmented in XY. Segmented cells are then measured using region_props. 
INPUTS:
o df: pandas DataFrame – A DataFrame with 5 columns and k rows where k is the number of segmented objects from the associated image. The columns from left to right are nuclear centroid position X, Y, Z and segmented nuclear volume. The fifth column should only have 6 rows: [maximum x coordinate, maximum y coordinate, maximum z coordinate, minimum x coordinate, minimum y coordinate, minimum z coordinate].
o image: 4D numpy array – An array with shape ~mxnx512x512 where m is the number of slices (minimum 5) and n is the number of channels (maximum 4) which can be in any order.  
o plot: boolean, optional – Default option is to not plot the nuclear distribution, but will if True or if save_name is defined.   
o save_name: string, optional – A string specifying the name to save a dictionary of individual cell properties if desired (labels, areas, centroid rows, centroid columns, perimeters, eccentricities, circularities) as a .csv file. Defaults to not saving this dictionary.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o average_cell_area: float – Average cell area in µm2.
o average_cell_perimeter: float – Average cell area in µm.
o average_cell_circularity: float – Average cell circularity.
 
#### batch_process (list_of_dfs, list_of_unshuffled_images, list_of_names, save_name = False, actin_channel = 1, invert = False)
This function allows for processing multiple images at once. Detailed descriptions of each image will be output into a pandas DataFrame, and are optionally saved to a .csv file. 
INPUTS:
o list_of_dfs: list of pandas DataFrames – all DataFrames associated with their respective images that will be analyzed  
o list_of_unshuffled_images: list of 4D numpy arrays – all images that will be analyzed
o list_of_names: list of strings – the names of all images that will be analyzed  
o save_name: string, optional – A string specifying the name to save a dictionary of analyzed image properties, defined below as a .csv file. Defaults to not saving this dictionary.
o actin_channel: int, optional – An integer to specify which of the channels for a multichannel image corresponds to the actin channel. Default is set to 1, the expected DAPI channel is 0. 
o invert: boolean, optional – Option to invert the z slices of an image if the z-stack was taken top to bottom instead of bottom to top. Defaults to False for images imaged bottom to top.   
OUTPUTS:
o analyzed_df: pandas DataFrame – DataFrame consisting of columns: image_name, layer_classifications, cells_in_layer, cells_above_layer, total_cells, cell_density, layer_height, percent_above, average_cell_area, average_cell_perimeter, average_cell_circularity, equivalences (Lateral-to-Apical shape index), derivative_peak_rule (True if No Shoulder). 

This readme was written by Drs Christian Cammarota and Tara M Finegan in September 2023.
