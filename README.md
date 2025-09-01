This folder contains the code used in the project “Sample Efficient Likelihood-Free Inference for Virus Dynamics with Different Types of Experiments” by Yingying Xu and Ulpu Remes.
Paper: https://doi.org/10.48550/arXiv.2508.11042




File Explanation

Data folder: ‘psimon-sH1N1’ (containing four datasets) 
The original source is Simon, P. F. et al., (2016), titled “Avian influenza viruses that cause highly virulent infections in humans exhibit distinct replicative properties in contrast to human H1N1 viruses”, Scientific Reports, Vol. 6, p. 24154. doi: 10.1038/srep24154. Available at: https://www.nature.com/articles/srep24154.

Simulator components: ‘simulator.py’ and ‘meanfield_solve_ivp_w_vrna.py’  (Note: ‘simulator.py’ requires ‘meanfield_solve_ivp_w_vrna.py’ to run) 
The simulator in this project combines components from the original sources by Quirouette, Christian et al., (2023). “The effect of random virus failure following cell entry on infection outcome and the success of antiviral therapy”, Scientific Reports, Vol. 13 No. 1, p. 17243. doi: 10.1038/s41598-023-44180-w.

Inference using real experimental data: ‘MFM_BOLFI_fail_detection_multi_Nov14.ipynb’ 

Inference using synthetic data: ‘MFM_BOLFI_fail_detection_multi_synthetic_2_Aug25.ipynb'

Plotting inference results for real experimental data as a corner plot: ‘plotresult.py’

Plotting inference results for synthetic data: ‘plot_synthetic.ipynb'

Saved learned Gaussian process model for inference using real experimental data: virology_GPmodel_bolfi_n1000_Nov14.npz  (corresponds to the ‘MFM_BOLFI_fail_detection_multi_Nov14.ipynb’ file)

Saved learned Gaussian process model for inference using synthetic data: virology_GPmodel_bolfi_synthetic_n1000_wider_Aug25.npz (corresponds to the  ‘MFM_BOLFI_fail_detection_multi_synthetic_2_Aug25.ipynb’ file)

Saved samples from posterior distribution in inference using real experimental data: sample_nsample100000_result_bolfi_fail_detecton_multi_euclidean_data_n1000_Nov14.pkl (The   raw file is large; a compressed zip file is uploaded)



A demo notebook with a 2D toy example demonstrating how BOLFI with classifier algorithm works:
Remes, Ulpu (2024). “ELFI Notebooks: Failure-robust BOLFI”, Accessed: 2024-2-20. available at:
https://github.com/uremes/elfi-notebooks/blob/extend_model/failure_robust_BOLFI.ipynb.
