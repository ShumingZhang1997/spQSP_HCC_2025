#include "DiffuseGrid.h"

#include "../core/GlobalUtilities.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

DiffuseGrid::DiffuseGrid()
{
	setup_biofvm_grid();
}


DiffuseGrid::~DiffuseGrid()
{
}

/*! initializae tumor microenvironment grid.
    This is called during model setup.
	No need to serialize variables that remains unchanged 
	during simulation.
*/
void DiffuseGrid::setup_biofvm_grid(void) {

	_tme.name="tumor_pde";

	// density unit: amount/cm^3. 'amount' correspond to 
	// amount/sec in release rate unit.
	_tme.set_density(0, "IFNg" , "amount/mL" );
	_tme.add_density("IL_2", "amount/mL");
	_tme.add_density("CCL2", "amount/mL");
	_tme.add_density("ArgI", "amount/mL");
	_tme.add_density("NO", "amount/mL");	
	_tme.add_density("TGFB", "amount/mL");
	_tme.add_density("IL10", "amount/mL");
	_tme.add_density("IL12", "amount/mL");
	_tme.add_density("VEGFA", "amount/mL");
	_tme.add_density("O2", "amount/mL"); // O2 using unit of mg/mL
	//_tme.add_density("IL_10", "amount/micron^3");

	_tme.spatial_units = "cm";
	_tme.mesh.units = "cm";
	_tme.time_units = "sec";

	double minX, minY, minZ, maxX, maxY, maxZ, mesh_resolution;

	minX = 0;
	minY = 0;
	minZ = 0;
	// convert from micron to cm
	mesh_resolution = params.getVal(PARAM_VOXEL_SIZE_CM);
	maxX = params.getVal(PARAM_TUMOR_X) * mesh_resolution;
	maxY = params.getVal(PARAM_TUMOR_Y) * mesh_resolution;
	maxZ = params.getVal(PARAM_TUMOR_Z) * mesh_resolution;

	_tme.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	//_tme.display_information( std::cout );

	// register the diffusion solver 	
	_tme.diffusion_decay_solver = BioFVM::diffusion_decay_solver__constant_coefficients_LOD_3D; 
	
	// register substrates properties 
	_tme.diffusion_coefficients[CHEM_IFN] = params.getVal(PARAM_IFN_G_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_IFN] =  params.getVal(PARAM_IFN_G_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_IL_2] = params.getVal(PARAM_IL_2_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_IL_2] =  params.getVal(PARAM_IL_2_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_CCL2] = params.getVal(PARAM_CCL2_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_CCL2] =  params.getVal(PARAM_CCL2_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_ARGI] = params.getVal(PARAM_ARGI_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_ARGI] =  params.getVal(PARAM_ARGI_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_NO] = params.getVal(PARAM_NO_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_NO] =  params.getVal(PARAM_NO_DECAY_RATE); // 1/sec	
	
	_tme.diffusion_coefficients[CHEM_TGFB] = params.getVal(PARAM_TGFB_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_TGFB] = params.getVal(PARAM_TGFB_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_IL_10] = params.getVal(PARAM_IL_10_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_IL_10] =  params.getVal(PARAM_IL_10_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_IL_12] = params.getVal(PARAM_IL_12_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_IL_12] = params.getVal(PARAM_IL_12_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_VEGFA] = params.getVal(PARAM_VEGFA_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_VEGFA] = params.getVal(PARAM_VEGFA_DECAY_RATE); // 1/sec

	_tme.diffusion_coefficients[CHEM_O2] = params.getVal(PARAM_O2_DIFFUSIVITY); // cm^2 / sec 
	_tme.decay_rates[CHEM_O2] = params.getVal(PARAM_O2_DECAY_RATE); // 1/sec


	_nrVoxel = _tme.number_of_voxels();
	_nrSubstrate = _tme.number_of_densities();

	for( int i=0; i < _tme.number_of_voxels() ; i++ )
	{
		_tme.density_vector(i)[0]= 0.0;
		_tme.density_vector(i)[CHEM_O2] = 0.673;
	}
	for (size_t i = 0; i < _nrSubstrate; i++) {
		_min_substrate.push_back(1E-5);
	}

	return;
}
/*! check if tme can be skipped
    Skip time step If:
	1. total amount of substrate is smaller than set value
	2. No existing source/sink
*/
bool DiffuseGrid::grid_skippable(void) {
	bool res = true;
	for (size_t i = 0; i < _nrSubstrate; i++)
	{
		res &= (get_total_amount(i) < _min_substrate[i]);
	}
	res &= _sink_source.empty();
	return res;
}

};
};
