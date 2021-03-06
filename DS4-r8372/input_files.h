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
	
	input_files.push_back(prefix + "8372.0-99_ntuple.root");
	input_files.push_back(prefix + "8372.100-199_ntuple.root");
	input_files.push_back(prefix + "8372.200-299_ntuple.root");
	input_files.push_back(prefix + "8372.300-399_ntuple.root");
	input_files.push_back(prefix + "8372.400-499_ntuple.root");
	input_files.push_back(prefix + "8372.500-599_ntuple.root");
	input_files.push_back(prefix + "8372.600-699_ntuple.root");
}
