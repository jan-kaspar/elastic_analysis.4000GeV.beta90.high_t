#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

std::string input_ntuple_name;

void InitInputFiles()
{
	input_ntuple_name = "TotemNtuple";

	input_files.clear();

	// 2012_07_12-13_b_hsx

	std::string prefix = "rfio:///castor/cern.ch/totem/offline/Reco/2012/Physics/July/v_3.11/";
	
	/*
	input_files.push_back(prefix + "8369.0-99_ntuple.root");
	input_files.push_back(prefix + "8369.100-199_ntuple.root");
	input_files.push_back(prefix + "8369.200-299_ntuple.root");
	input_files.push_back(prefix + "8369.300-399_ntuple.root");
	input_files.push_back(prefix + "8369.400-499_ntuple.root");
	input_files.push_back(prefix + "8371.0-99_ntuple.root");
	input_files.push_back(prefix + "8371.100-199_ntuple.root");
	input_files.push_back(prefix + "8371.200-299_ntuple.root");
	input_files.push_back(prefix + "8371.300-399_ntuple.root");
	*/
	input_files.push_back(prefix + "8372.0-99_ntuple.root");
	input_files.push_back(prefix + "8372.100-199_ntuple.root");
	input_files.push_back(prefix + "8372.200-299_ntuple.root");
	input_files.push_back(prefix + "8372.300-399_ntuple.root");
	input_files.push_back(prefix + "8372.400-499_ntuple.root");
	input_files.push_back(prefix + "8372.500-599_ntuple.root");
	input_files.push_back(prefix + "8372.600-699_ntuple.root");

	/*
	std::string prefix = "rfio:///castor/cern.ch/totem/offline/Reco/2012/Physics/July/v_3.8/validation_split/";
	
	input_files.push_back(prefix + "val8369aa/val8369aa_totem_ntuple.root");
	input_files.push_back(prefix + "val8369ab/val8369ab_totem_ntuple.root");
	input_files.push_back(prefix + "val8369ac/val8369ac_totem_ntuple.root");
	input_files.push_back(prefix + "val8369ad/val8369ad_totem_ntuple.root");
	input_files.push_back(prefix + "val8369ba/val8369ba_totem_ntuple.root");
	input_files.push_back(prefix + "val8369bb/val8369bb_totem_ntuple.root");
	input_files.push_back(prefix + "val8369bc/val8369bc_totem_ntuple.root");
	input_files.push_back(prefix + "val8369bd/val8369bd_totem_ntuple.root");
	input_files.push_back(prefix + "val8369ca/val8369ca_totem_ntuple.root");
	input_files.push_back(prefix + "val8369cb/val8369cb_totem_ntuple.root");
	input_files.push_back(prefix + "val8369cc/val8369cc_totem_ntuple.root");
	input_files.push_back(prefix + "val8369cd/val8369cd_totem_ntuple.root");
	input_files.push_back(prefix + "val8369da/val8369da_totem_ntuple.root");
	input_files.push_back(prefix + "val8369db/val8369db_totem_ntuple.root");
	input_files.push_back(prefix + "val8369dc/val8369dc_totem_ntuple.root");
	input_files.push_back(prefix + "val8369dd/val8369dd_totem_ntuple.root");
	input_files.push_back(prefix + "val8369ea/val8369ea_totem_ntuple.root");
	input_files.push_back(prefix + "val8369eb/val8369eb_totem_ntuple.root");
	input_files.push_back(prefix + "val8369ec/val8369ec_totem_ntuple.root");
	
	input_files.push_back(prefix + "val8371aa/val8371aa_totem_ntuple.root");
	input_files.push_back(prefix + "val8371ab/val8371ab_totem_ntuple.root");
	input_files.push_back(prefix + "val8371ac/val8371ac_totem_ntuple.root");
	input_files.push_back(prefix + "val8371ad/val8371ad_totem_ntuple.root");
	input_files.push_back(prefix + "val8371ba/val8371ba_totem_ntuple.root");
	input_files.push_back(prefix + "val8371bb/val8371bb_totem_ntuple.root");
	input_files.push_back(prefix + "val8371bc/val8371bc_totem_ntuple.root");
	input_files.push_back(prefix + "val8371bd/val8371bd_totem_ntuple.root");
	input_files.push_back(prefix + "val8371ca/val8371ca_totem_ntuple.root");
	input_files.push_back(prefix + "val8371cb/val8371cb_totem_ntuple.root");
	input_files.push_back(prefix + "val8371cc/val8371cc_totem_ntuple.root");
	input_files.push_back(prefix + "val8371cd/val8371cd_totem_ntuple.root");
	input_files.push_back(prefix + "val8371da/val8371da_totem_ntuple.root");
	input_files.push_back(prefix + "val8371db/val8371db_totem_ntuple.root");
	input_files.push_back(prefix + "val8371dc/val8371dc_totem_ntuple.root");
	input_files.push_back(prefix + "val8371dd/val8371dd_totem_ntuple.root");
	
	input_files.push_back(prefix + "val8372aa/val8372aa_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ab/val8372ab_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ac/val8372ac_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ad/val8372ad_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ba/val8372ba_totem_ntuple.root");
	input_files.push_back(prefix + "val8372bb/val8372bb_totem_ntuple.root");
	input_files.push_back(prefix + "val8372bc/val8372bc_totem_ntuple.root");
	input_files.push_back(prefix + "val8372bd/val8372bd_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ca/val8372ca_totem_ntuple.root");
	input_files.push_back(prefix + "val8372cb/val8372cb_totem_ntuple.root");
	input_files.push_back(prefix + "val8372cc/val8372cc_totem_ntuple.root");
	input_files.push_back(prefix + "val8372cd/val8372cd_totem_ntuple.root");
	input_files.push_back(prefix + "val8372da/val8372da_totem_ntuple.root");
	input_files.push_back(prefix + "val8372db/val8372db_totem_ntuple.root");
	input_files.push_back(prefix + "val8372dc/val8372dc_totem_ntuple.root");
	input_files.push_back(prefix + "val8372dd/val8372dd_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ea/val8372ea_totem_ntuple.root");
	input_files.push_back(prefix + "val8372eb/val8372eb_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ec/val8372ec_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ed/val8372ed_totem_ntuple.root");
	input_files.push_back(prefix + "val8372fa/val8372fa_totem_ntuple.root");
	input_files.push_back(prefix + "val8372fb/val8372fb_totem_ntuple.root");
	input_files.push_back(prefix + "val8372fc/val8372fc_totem_ntuple.root");
	input_files.push_back(prefix + "val8372fd/val8372fd_totem_ntuple.root");
	input_files.push_back(prefix + "val8372ga/val8372ga_totem_ntuple.root");
	input_files.push_back(prefix + "val8372gb/val8372gb_totem_ntuple.root");
	input_files.push_back(prefix + "val8372gc/val8372gc_totem_ntuple.root");
	input_files.push_back(prefix + "val8372gd/val8372gd_totem_ntuple.root");
	*/
}
