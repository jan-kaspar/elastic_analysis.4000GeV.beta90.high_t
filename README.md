dataset: DS4
final output: file=DS-merged/merged.root, object=bt1/DS4-sc/combined/h_dsdt TODO: confirm

--------------------

1) distill: input=original ntuples, output=distill_<dgn> (distilled ntuples)

2) distributions: input=distill_<dgn>, output=distributions_<dgn>
	initial round
	use_time_dependent_resolutions = false
	use_3outof4_efficiency_fits = false 

3) resolution fit: input=distributions_<dgn>, output=resolution_fit_<dgn>

4) link eff3outof4_details_fit.root from main analysis as eff3outof4_details_fit_old.root

5) distributions: input=distill_<dgn>, output=distributions_<dgn>
	use_time_dependent_resolutions = true 
	use_3outof4_efficiency_fits = true

--------------------
unfolding

12) combine_distributions: input=distributions_<dgn>, output=combine_distributions

13) unfolding_cf: input=distributions_<dgn>, combine_distributions, smearing from anal, output=unfolding_fr

14) smearing_matrix_mc: input=smearing from anal, output=smearing_matrix_<dgn>

15) unfolding_gr: input=distributions_<dgn>, smearing_matrix_mc_<dgn>, output=unfolding_gr

16) re-run distributions to apply unfolding correction from 13)

17) DS-merged/merge: merge results from 17)

--------------------
systematics

20) models/sample_models: input=DS4/unfolding_cf_<dgn> output=model

21) systematics_ni: input=model (point 20), output=systematics_ni_<dgn>

22) systematics_matrix_mc: input=systematics_matrix_<dgn>, output=systematics_matrix




