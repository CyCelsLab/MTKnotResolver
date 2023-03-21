%% Written By Dhruv Khatri, IISER Pune, Cycels Lab  
clear; close all; 
inputFile = "61 - dyn1-20002 (copy).nd2 (series 1).mat"; 
A = load(inputFile) % Load the data_struct 
data_struct = A.data_struct
resolveCoordinates = SingleFilamentResolver(data_struct, frameIgnore, frameEdit, resfullPath); 
