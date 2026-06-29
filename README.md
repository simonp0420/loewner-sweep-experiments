# Original Description
Implementation of adaptive and semi-adaptive sweep algorithms for broadband electromagnetic simulation using the Loewner framework.

Based on:

T. N. Shilpa and R. Sinha, "Fully Adaptive and Semi-Adaptive Frequency Sweep Algorithm Exploiting Loewner-State Model for EM Simulation of Multiport Systems," IEEE Transactions on Microwave Theory and Techniques, vol. 73, no. 9, pp. 6245-6259, Sept. 2025, doi: 10.1109/TMTT.2025.3557208.

S. T. N and R. Sinha, "An Adaptive Frequency Sweep Algorithm for Broadband EM Simulation Incorporating Random Matrix Transformation in the Loewner Framework," 2025 IEEE Microwaves, Antennas, and Propagation Conference (MAPCON), Kochi, India, 2025, pp. 1-4, doi: 10.1109/MAPCON65020.2025.11426247.

The code may later be registered as a standalone module, or included in a sweep library with several methods, if I have inspiration and time.

The directory [examples_pssfss](examples_pssfss) contains examples using [PSSFSS.jl](https://github.com/simonp0420/PSSFSS.jl). Suggestions for other EM simulation libraries where these sweep methods could be tested are welcome.

<img width="1120" height="620" alt="image" src="https://github.com/user-attachments/assets/9ff2a661-3295-4ec8-b3a3-8b4b3ac542d5" />


# Updates from psimon_benchmarks Branch
Thanks to this amazing testing framework created by @uvegege it was easy to run a number of comparisons, which are described here.  Here is a summary of the 
timing and accuracy results obtained over all 10 test cases using the original default settings (shown near the bottom of the plot):

<img width="1000" height="620" alt="FirstRun" src="https://github.com/user-attachments/assets/1f630af8-7e2d-44a5-a916-8baac522deac" />

The error plot in the top half of the figure shows that three of the cases exhibited accuracy judged to be unacceptable (greater than 1e-4) for both Loewner and Random runs.
In the timing plot (bottom half), the gray dashed line
at ordinate value 1 represents the time required for the PSSFSS fast sweep.  Points plotted above this line indicate cases that were slower than the PSSFSS fast
sweep.  The Jerusalem Cross case using the Loewner algorithm was about 50 times slower!  This is because of the extensive linear algebra calculations needed
and the use of OpenBLAS on myy Intel machine.  The same case run using the MKL linear algebra library produced these results:

<img width="1000" height="620" alt="MKL" src="https://github.com/user-attachments/assets/535e51f9-2f12-4691-ac78-b20a16fde36b" />
 
Here we see that the execution times are more reasonable. For most cases, Loewner and Random are faster than the reference (PSSFSS Fast Sweep), but they are 
still both slower for the resistive square patch.  The accuracy for Loewner and Random are judged to be inadequate for Flexible Absorber, Jerusalem Cross,
and Resistive Square Patch.

Next is a case where the tolerance parameter `tol` was reduced to 0.0001 from its previous value of 0.001, and the number of parallel tasks `nthreads` was reduced to 
14 from its previous value of 24.  This latter change was found to have no effect and was retained in all the following cases reported here.

<img width="1000" height="620" alt="MKL_14tasks_tol1e-4" src="https://github.com/user-attachments/assets/44678ab2-ce0e-4694-98f1-4554d735e869" />

The accuracy is improved for some cases, but not for the ones previously judged inadequate.

Because the authors of the cited references used nonsingular (full-rank) random matrices in their Random formulation, I implemented this change (see src/RandomLU.jl)
and reverted to using LU decomposition rather than QR in hopes of obtaining better results in terms of Random accuracy and timing.  Here is a case with
tol = 1e-4 and nonsingular random matrices:

  <img width="1000" height="620" alt="MKL_14tasks_tol1e-4_complexWaWbLU" src="https://github.com/user-attachments/assets/877d38a0-15a7-4293-9237-306c625d281c" />

  Neither the speed nor accuracy of the Random case seemed to be very strongly affected by this change.

  The following results use the "semi-adaptive" formulation with `tol = 1e-3`, `l = 0.02` and `p = 4`:

  <img width="1000" height="620" alt="MKL_semi_14tasks_tol1e-3_complexWaWbLU" src="https://github.com/user-attachments/assets/ebf6cdae-ce3a-4494-87b7-484c0a61429c" />

The errors for the Loewner runs are borderline acceptable, though 2 of the Loewner times are excessive.

Next, we double the `l` parameter to `0.04` to force more initial samples, hoping for even greater accuracy:

<img width="1000" height="620" alt="MKL_semi_14tasks_tol1e-3_l0 04_complexWaWbLU" src="https://github.com/user-attachments/assets/d15b5707-9e51-44a7-a309-121d9d37d8a7" />

Here most cases show acceptably small errors, with some exceptions.  For most of the cases, the new algorithms are modestly faster than PSSFSS Fast Sweep, but both are
extremely slow for Flexible Absorber.

# Conclusion
I haven't been able to choose settings such that the new algorithms are reliably accurate and faster than (or at least not slower than) the presently implemented PSSFSS Fast Sweep algorithm.  So they will not be implemented at this time in PSSFSS.


  
