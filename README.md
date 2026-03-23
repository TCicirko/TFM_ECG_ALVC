# TFM_ECG_ALVC
ECG signal processing code
This MATLAB script performs automated longitudinal ECG analysis for patients with Arrhythmogenic Left Ventricular Cardiomyopathy (ALVC) using ECGkit.
The workflow is divided into three main parts.
Part 1 processes each ECG recording individually. It loads the filtered 12-lead ECG, detects R peaks using the Pan–Tompkins algorithm, extracts median beats for each lead, builds a synthetic 10-beat signal, performs waveform delineation, attempts recovery of missing P-wave markers from the previous beat when possible, and computes lead-wise biomarkers. The main biomarkers saved are QRS duration, PR interval, QT interval, QRS peak-to-peak amplitude, maximum and minimum QRS peak, QRS energy, average QRS power, and absolute QRS area. Outputs are saved in the patient/date folders, mainly inside QRSDetection, Delineation, and Markers.
Part 2 selects patients for longitudinal comparison. For each patient, the earliest usable year is treated as Year 1. The script then searches for the first usable recording in baseline year + 2 for the Year 1 vs Year 3 comparison, and baseline year + 4 for the Year 1 vs Year 5 comparison. A recording is considered usable if core biomarkers are available in at least 8 leads. Selection logs are saved as Excel files in the main PatientSamples folder.
Part 3 performs the statistical analysis and plotting. For each included patient, the script compares Year 1 with the follow-up recording, generates lead-wise bar plots and box plots for each biomarker, and stores them in patient-specific result folders. For group statistics, the median biomarker value across valid leads is used for each patient. Normality of patient-level differences is checked using the Lilliefors test. If any biomarker fails normality, all biomarkers in that comparison are tested with the Wilcoxon signed-rank test; otherwise, paired t-tests are used. Final statistical results are saved as Excel files in the PatientSamples folder.
Required folder structure
The script expects:
one folder per patient inside PatientSamples
one subfolder per ECG date, named as yyyymmdd
inside each date folder:
patient_dateBRLPF.mat
patient_date_HEADER.mat
Main output files
The script generates:
QRS detection plots
median beat plots
beat-3 marker plots
biomarker .mat files
P-marker rescue log
Year 1 vs Year 3 selection log
Year 1 vs Year 5 selection log
statistical result spreadsheets
patient-specific comparison plots
Important notes
QT is computed during processing but excluded from the final longitudinal statistical comparisons. Biomarkers are only computed for leads where the required markers are detected in physiologically valid positions. The script is designed for reproducible automated analysis, but results should still be reviewed visually, especially for cases with rescued P-wave markers or difficult delineation.