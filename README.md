# RHEP-JOC
Data and Code used to generate numerical results in the paper "A Computational Approach to First Passage Problems of Reflected Hyper-Exponential Jump Diffusion Processes," which is forthcoming at the INFORMS Journal on Computing.

All algorithms are coded in MATLAB 2017a. All data and intermediate results are store in .mat format of MATLAB. 

Users should first install the MATLAB "Multiprecision Computing Toolbox," which is available at https://www.advanpix.com.

To replicate the numerical results reported in Section 4 of the paper, please use "Figure1_Q2_to_Q5.m" to generate Figure 1; use "Table3_PanelA_Dist_FPT.m" to generate LIVs in Panel A of Table 3; use "Table3_PanelB_Joint_Dist_FPT_TV.m" to generate LIVs in Panel B of Table 3; use "Table5_Table6_Q5.m" to generate Table 5 (set mydigits=200) and Table 6 (set mydigits=32). To generate Monte Carlo values, one can use "RHEP_path.m" to generate sample paths of RHEPs. 

For any question about the use of these data and code, please contact Xuewei Yang, the corresponding author, via xwyang@nju.edu.cn or xwyang@aliyun.com
