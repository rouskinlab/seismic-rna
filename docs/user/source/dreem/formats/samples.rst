Samples
+++++++


**Example: samples.csv**

 ================ ======== ======== ========== =============== =================== ============= ======== =========== 
  sample           user     date     exp_env    temperature_k   inc_time_tot_secs   DMS_conc_mM   buffer   cell_line  
 ================ ======== ======== ========== =============== =================== ============= ======== =========== 
  01_1_S22_reads   Lauren   1/5/23   in_vitro   310             300                 10.5          buffer              
 ================ ======== ======== ========== =============== =================== ============= ======== =========== 


**About:** Add information related to the entire sample set by writing a `samples.csv` file, which is a tab-delimited text file with the following columns:

===================================== ======================================================================================================================================= 
  column                                description                                                                                                                                                           
===================================== ======================================================================================================================================= 
  sample                                sample name, corresponding to the prefix of the fastq file (ex: sample 'sample1' corresponds to the input file 'sample1_R1.fastq')  
  user                                  user name                                                                                                                          
  date                                  date of the experiment                                                                                                             
  exp_env                               experimental environment, in_vivo or in_vitro                                                                                      
  temperature_k                         temperature in Kelvin                                                                                                              
  inc_time_tot_secs                     total incubation time in seconds                                                                                                   
  DMS_conc_mM                           concentration of DMS in mM                                                                                                         
  buffer (if exp_env is in_vitro)       exact buffer including Mg, eg 300mM Sodium Cacodylate , 3mM Mg                                                                     
  cell_line (if exp_env is in_vivo)     cell line                                                                                                                          
  [any other columns]                   any other columns will be added to corresponding sample in the output file                                                             
===================================== ======================================================================================================================================= 
