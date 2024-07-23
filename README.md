# KnotResolver: MATLAB Code Guide

![A beating filament from a pinned gliding assay](Images/density27len47(MOD_trans)contour.gif)

KnotResolver is a MATLAB script that can be used to segment and track single filaments. The script is designed for a “beating assay” setup.

## Table of Contents
- [Download source code](#download-source-code)
- [Prerequisites](#prerequisites)
- [1. Workflow](#1-workflow)
  - [1.1. Creating a ‘.csv’ file for analysis](#11-creating-a-csv-file-for-analysis)
  - [1.2. Segmentation parameters](#12-segmentation-parameters)
  - [1.3. KnotResolver.m](#13-knotresolverm)
  - [1.4. Output](#14-output)
    - [DemoOutput](#demooutput)

## Download source code
The source code can be downloaded from GitHub: [https://github.com/CyCelsLab/MTKnotResolver](https://github.com/CyCelsLab/MTKnotResolver)

## Prerequisites
In order to use KnotResolver via source code, the user needs to install MATLAB with the following toolboxes:

- MATLAB 9.12
- Signal Processing Toolbox 9.0
- Image Processing Toolbox 11.5
- Statistics and Machine Learning Toolbox 12.3
- Curve Fitting Toolbox 3.7
- Parallel Computing Toolbox 7.6 (Optional)
- MATLAB Parallel Server 7.6 (Optional)
- Polyspace Bug Finder 3.6 (Optional)

The interface of KnotResolver requires MATLAB version R2019b or higher.

Note: Please close other figure tabs if open in MATLAB, before running KnotResolver.

## 1. Workflow

### 1.1. Creating a ‘.csv’ file for analysis
The analysis begins by collating all time series containing a single microtubule pinned at one end. The details required in the `.csv` file headers are "FileName" and "FolderPath." These headers are case-sensitive and should reflect the name of the time series with extensions and the location of the parent folder. The code expects input image time series in 8-bit `.tif` format. An example `.csv` file and two default time series are available on GitHub.

Since different images can contain differences in noise and foreground, it is important to tune the parameters of segmentation. This can be done interactively using the `KnotSegmenter.m` file. Simply run the MATLAB file and tune the parameters until satisfied with the segmentation results.

The output is a `.csv` file with the suffix "optimized." This file contains the optimized segmentation parameters that can be run using the `KnotResolver.m` file. This `.csv` will be saved in the same directory as the original `.csv` file.

### 1.2. Segmentation parameters
- **Threshold**: Threshold tuning between 0 and 1. Used as an input for the `imbinarize` intensity threshold method.
- **Contraction Bias and Smooth Factor**: Input parameters for the `activecontour` method. Contraction bias ranges from -1 to +1. Positive values force the contours to shrink, and vice versa. The smooth factor determines the smoothness of the optimized curve.
- **Iterations**: Number of iterations for the `activecontour` method.

### 1.3. KnotResolver.m
The `.csv` with the optimized parameters can be run directly with the `KnotResolver.m` script. For each file, the code shows the segmentation overlay on the original input image with a crosshair for MT selection. The central microtubule needs to be selected. The contour selection location determines the position of the fixed tip. Therefore, careful selection of the contour near the fixed tip should be made.

Once the selection of the central contour is made, the MT in the next frame is identified based on overlap with the initial selection, which is updated in every frame.

### 1.4. Output

#### DemoOutput
A `DemoOutput` is created if the `InputFile_optimized.csv` file is used to run the `KnotResolver.m` script. The folder contains the segmentation output in a `.mat` file which has the same name as the input image. This `.mat` file contains the information of image size and the binary segmented skeleton. This file is created using the optimized segmentation parameters.

A separate folder with the same name as the input image is created and contains the following outputs:

- **ContourOverlay.pdf**: Contains a subplot of three figures arranged vertically. The first figure contains the overlay of the complete resolved time series color coded along the length. 
- **TangentAngleKymograph.pdf**: Contains the kymograph of the tangent angles along the microtubule contours.
- **End2End.pdf**: Contains the distance between the two ends of the filament in the time series.
- **MTdata.csv**: Contains three columns that correspond to the frame number, X and Y coordinate of the skeletons. The indices are in pixels and need to be rescaled with the correct scaling factor.
- **resolveCoordinates_v*.mat**: `.mat` file containing the resolved coordinates. This can be used for further analysis when required. The contour indices are arranged from tip to bottom along with the frame number.
- **segOverlay_v*.tif**: Final output skeleton color coded along the length. This can be visually checked to confirm if the output is correct.
- **TipAngle.pdf**: FFT analysis of the tip angle suspended at the two tips of the beating microtubule. The subplots show the tip angle change with frame number and power spectrum output from the FFT analysis.
- **statsSummary.txt**: This file contains the log of input parameters of scaling, time step, and ignored distance from the tip. Additionally, dominant frequency and contour length parameters with frame-wise change are also documented.
