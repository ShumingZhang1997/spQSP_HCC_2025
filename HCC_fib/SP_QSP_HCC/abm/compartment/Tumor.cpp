//#include <boost/serialization/export.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
#include "Tumor.h"
#include "LymphCentral.h"
//BOOST_CLASS_EXPORT_IMPLEMENT(Tumor);

//#include "../agent/CancerCell.h"`
#include <algorithm>
#include <exception>
#include <numeric>
#include <queue>
#include <cmath>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include "../../core/GlobalUtilities.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//#include "../../MemoryUsage.h"


Tumor::Tumor(int x, int y, int z)
	: SpatialCompartment(x, y, z)
	, _stats()
	, _voxel_ic()
	, _chem()
	, _init_vas()
	, _t_source()
	, _mdsc_source()
	, _mac_source()
	, _tInitDummy(new TCell(this))
	, _cInitDummy(new CancerCell(this))
	, _macInitDummy(new Mac(this))
	, _fibInitDummy(new Fib(this))
	, _tcd4InitDummy(new TCD4(this))
	, _MDSCInitDummy(new MDSC(this))
	, _VasInitDummy(new Vas(this))
	//, _voxel_size(params.getVal(PARAM_VOXEL_SIZE))
	, _allow_shift_grid(false)
	, _center_target()
	, _agGrid_temp(AgentGrid(x, y, z, NULL))
	, _var_abm_to_qsp()
	, _concentration_cc(0)
	, _concentration_t_cyt(0)
	, _concentration_t_reg(0)
	, _concentration_t_h(0)
	, _concentration_t_eff_tum(0)
	, _concentration_m1(0)
	, _concentration_m2(0)
	, _concentration_mdsc(0)
	, _concentration_nivo(0)
	, _concentration_ipi(0)
	, _concentration_ent(0)
	, _concentration_cabo(0)
	, _concentration_cx(0)
	, _concentration_t_exh(0)
	, _tumor_volume(0)
	, _f_tum_cap(0)
	, _R_cabo(0)
	, _cancer_counts(0)
	, _next_unique_id(0)
{
	//CC total, CC death total, CC death Teff, Teff recruit, TCD4 recruit
	_var_abm_to_qsp = std::vector<double>(TUMEX_VAR_NUM, 0);

	// dummy cells do not act as source
	//_tInitDummy->get_source_IFNg()->set_remove();
	initAgentGrid();
}

/*!
  _agGrid uses pointer to AgentGridVoxel (or its derived classes)
  as its element. This has a higher computational cost than directly
  hosting the instantications of the derived class (base won't work).
  The benefit is that this way the voxel objects can be accesses from
  the SpatialCompartment abstract class.
*/
bool Tumor::initAgentGrid() {
	//std::cout << "init agGrid:" << std::endl;
	for (int i = 0; i < _sizeX; i++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int k = 0; k < _sizeZ; k++)
			{
				_agGrid(i, j, k) = new TumorGridVoxel(i, j, k);
				//std::cout << _agGrid(i, j, k) << std::endl;
			}
		}
	}

	return true;
}

Tumor::~Tumor()
{
	//std::cout << "destructor:" << std::endl;
	for (int i = 0; i < _sizeX; i++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int k = 0; k < _sizeZ; k++)
			{
				//std::cout << _agGrid(i, j, k) << std::endl;
				delete _agGrid(i, j, k);
			}
		}
	}

	delete _tInitDummy;
	delete _cInitDummy;
	delete _macInitDummy;
	delete _fibInitDummy;
	delete _tcd4InitDummy;
	delete _MDSCInitDummy;
	delete _VasInitDummy;
}

// add initial vasculature cell
void Tumor::add_init_vas(const Coord3D& c) {
	_init_vas.push_back(c);
}

void Tumor::update_ECM() {
	Grid3D<double> fibroblast_field(_sizeX, _sizeY, _sizeZ, 0);
    const int field_radius = 10;
    const double standard_deviation = 3;
    const int kernel_size = 2 * field_radius + 1;

    // Precompute Gaussian kernel
    std::vector<std::vector<std::vector<double>>> gaussian_kernel(kernel_size,
        std::vector<std::vector<double>>(kernel_size, std::vector<double>(kernel_size, 0.0)));
    const double variance = standard_deviation * standard_deviation;
    //const double normalizer = 1.0 / (std::pow(standard_deviation, 3) * std::sqrt(2 * 3.1415926));
	const double normalizer = 1.0 / (std::pow(standard_deviation, 3) * std::sqrt(2 * 3.1415926));
    // Debug prints
    for (int x = -field_radius; x <= field_radius; x++) {
        for (int y = -field_radius; y <= field_radius; y++) {
            for (int z = -field_radius; z <= field_radius; z++) {
                double dist_sq = x * x + y * y + z * z;  // Could be up to 300 for radius 10
                double exp_term = -dist_sq / (2*variance); 
                
                //std::cout << "dist_sq: " << dist_sq << ", exp_term: " << exp_term << std::endl;
                
                gaussian_kernel[x + field_radius][y + field_radius][z + field_radius] =
                    normalizer * std::exp(exp_term);
            }
        }
    }

    // Iterate over fibroblast cells
    for (const auto& agent : _vecAgent) 
	{
		auto p = dynamic_cast<Cell_Tumor*>(agent);
		if (p->getType() != AgentTypeEnum::CELL_TYPE_FIB) {
			continue;
		}

		Coord3D cell_coord = p->getCoord();
		double scale = (p->getState() == AgentStateEnum::FIB_CAF) ? 1.0 : 0.5;

		// Define bounds
		int x_min = std::max(0, (int)(cell_coord.x - field_radius));
		int x_max = std::min(_sizeX, (int)(cell_coord.x + field_radius + 1));
		int y_min = std::max(0, (int)(cell_coord.y - field_radius));
		int y_max = std::min(_sizeY, (int)(cell_coord.y + field_radius + 1));
		int z_min = std::max(0, (int)(cell_coord.z - field_radius));
		int z_max = std::min(_sizeZ, (int)(cell_coord.z + field_radius + 1));

		// Apply the kernel to the region
		for (int x = x_min; x < x_max; x++) {
			for (int y = y_min; y < y_max; y++) {
				for (int z = z_min; z < z_max; z++) {
					if (fibroblast_field.inGrid(x, y, z)) {  // Ensure valid grid coordinates
						int kernel_x = x - (int)cell_coord.x + field_radius;
						int kernel_y = y - (int)cell_coord.y + field_radius;
						int kernel_z = z - (int)cell_coord.z + field_radius;

						if (kernel_x >= 0 && kernel_y >= 0 && kernel_z >= 0 &&
							kernel_x < kernel_size && kernel_y < kernel_size && kernel_z < kernel_size) {
							double old_field = fibroblast_field.get(x, y, z);
							double value = scale * gaussian_kernel[kernel_x][kernel_y][kernel_z];
							fibroblast_field.set(x, y, z, old_field + value);
							//std::cout << "kernel_x: " << kernel_x << ", kernel_y: " << kernel_y << ", kernel_z: " << kernel_z << ", value: " << value << std::endl;
						}
					}
				}
			}
		}
	}
	

	for_each_grid_coord(true, true, true, [&](Coord& c) {
		// first calculate the ECM decay from last time point, each voxel registers the concentration of ECM 
		auto p_agGrid = dynamic_cast<TumorGridVoxel*>(_agGrid(c));
		//get the fibroblast field and tgfb 
		double fib_field = fibroblast_field.get(c.x, c.y, c.z);
		double TGFB = get_chem(c, CHEM_TGFB);
		double H_CAF_TGFB = (TGFB / (TGFB + params.getVal(PARAM_FIB_CAF_EC50)));
		// get the ECM concentration and amount
		double old_ECM_concentration = p_agGrid->getECMDensity();
		double old_ECM_amount = old_ECM_concentration * std::pow(params.getVal(PARAM_VOXEL_SIZE_CM), 3); // 
		
		double saturation = old_ECM_concentration/params.getVal(PARAM_FIB_ECM_SATURATION);
		saturation = (saturation < 1) ? saturation : 1;
		double ECM_secretion = fib_field *  (1 + H_CAF_TGFB) * params.getVal(PARAM_FIB_ECM_RELEASE_CAF) / 3 * (1 - saturation);
		// [ECM]_n+1 = [ECM]_n + k * delta_t - d * [ECM]_n
		// k is the secretion and d is the degradation rate. 
		double current_ECM_amount = old_ECM_amount * std::exp(-params.getVal(PARAM_SEC_PER_TIME_SLICE) * params.getVal(PARAM_FIB_ECM_DECAY_RATE));
		double new_ECM_amount = current_ECM_amount + ECM_secretion;
		double new_ECM_concentration = new_ECM_amount / std::pow(params.getVal(PARAM_VOXEL_SIZE_CM), 3);


		/*
		std::cout << "fib_field: " << fib_field
			<< ", TGFB: " << TGFB
			//<< ", max_num_cell: " << max_num_cell
			<< ", H_CAF_TGFB: " << H_CAF_TGFB
			<< ", old_ECM_concentration: " << old_ECM_concentration
			<< ", old_ECM_amount: " << old_ECM_amount
			<< ", current_ECM_amount: " << current_ECM_amount
			<< ", saturation: " << saturation
			<< ", new_ECM_amount: " << new_ECM_amount
			<< ", new_ECM_concentration: " << new_ECM_concentration << std::endl;
		

		
		std::cout << "fib_field: " << fib_field
			<< ", TGFB: " << TGFB
			<< ", H_CAF_TGFB: " << H_CAF_TGFB
			<< ", old_ECM_concentration: " << old_ECM_concentration
			<< ", degration rate: " << params.getVal(PARAM_SEC_PER_TIME_SLICE) * params.getVal(PARAM_FIB_ECM_DECAY_RATE)
			<< ", secretion: " << ECM_secretion
			<< ", new_ECM_amount: " << new_ECM_amount
			<< ", new_ECM_concentration: " << new_ECM_concentration << std::endl;
		*/

		// ECM concentration never go below the baseline value
		new_ECM_concentration = (new_ECM_concentration > params.getVal(PARAM_FIB_ECM_BASELINE) ? new_ECM_concentration : params.getVal(PARAM_FIB_ECM_BASELINE));


		p_agGrid->setECMDensity(new_ECM_concentration);
		p_agGrid->setFibField(fib_field);

		});

}


// add vasculature entry points for all immune cells (each cell have one chance)
void Tumor::update_vas(){
	
	int neighborhood_radius = 3;
	int neighborhood_size = neighborhood_radius * 2 + 1;
	//Construct an agent grid with size X+2n-1, Y+2n-1, Z+2n-1
	//The grid is used to store the number of certain cell type in the moore neighborhood with radius n of coordinate C
	//access number of cell around Coord C by calling neighborhood_cancer_counts(c.x + neighborhood_radius, c.y + neighborhood_radius, c.z + neighborhood_radius)
	//Much faster compare to naive method
	Grid3D<int> neighborhood_cancer_counts(_sizeX + neighborhood_size - 1, _sizeY + neighborhood_size - 1, _sizeZ + neighborhood_size - 1, 0);
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		int cancer_count = 0;
		_agGrid(c)->countNumAgentByType(AgentTypeEnum::CELL_TYPE_CANCER, cancer_count, false);
		neighborhood_cancer_counts(c) = cancer_count;
	});
	neighborhood_cancer_counts.window_counts_inplace(neighborhood_size);


	//count all voxel around coord c
	Grid3D<int> neighborhood_voxel_counts(_sizeX + neighborhood_size - 1, _sizeY + neighborhood_size - 1, _sizeZ + neighborhood_size - 1, 0);
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		int count = 0;
		if (_agGrid.inGrid(c)) {
			neighborhood_voxel_counts(c) = 1;
		}
		});
	neighborhood_voxel_counts.window_counts_inplace(neighborhood_size);
	

	//update the cabozantinib resistance
	update_R_cabo();

	CellVec::iterator last = stable_partition(_vecAgent.begin(), _vecAgent.end(), [&](CellAgent* p) {
		auto pC = dynamic_cast<Cell_Tumor*>(p);
		return (pC->getType() == AgentTypeEnum::CELL_TYPE_VAS);
	});// sortLiveCellAddToStats);

	// Calculate the number of partitioned elements
	int partitioned_count = std::distance(_vecAgent.begin(), last);
	std::cout << "Number of elements of type CELL_TYPE_VAS: " << partitioned_count << std::endl;

	double max_cancer = _sizeX * _sizeY * _sizeZ;
	double tumor_scaler = std::sqrt(1e5 * max_cancer / (_concentration_cc * params.getVal(PARAM_AVOGADROS)));
	// 8e-5 / _tumor_volume as the scaling factor accounting for the volume of the tumor in the QSP model
	//double tumor_scaler = std::sqrt(8e-5 / _tumor_volume);
	std::cout << "tumor_scaler: " << tumor_scaler << std::endl;
	// Optional: calculate the number of remaining elements
	int remaining_count = std::distance(last, _vecAgent.end());
	std::cout << "Number of remaining elements: " << remaining_count << std::endl;

	for (CellVec::iterator it = _vecAgent.begin(); it < last; it++)
	{
		auto p = dynamic_cast<Cell_Tumor*>(*it);
		if (p->getType() == AgentTypeEnum::CELL_TYPE_VAS && p->getState() == AgentStateEnum::VAS_PHALANX)
		{
			//Adding O2 source to grid for all matured vasculature cells
			auto pVas = dynamic_cast<Vas*>(p);
			Coord3D c_vas = pVas->getCoord();
			// Oxygen transport mechanism adopted from Sharan et al.
			//calculate oxygen transport rate based on vasulature density and oxygen level
			double pi = 3.1415926;
			double O2 = get_chem(c_vas, CHEM_O2);
			// dimensionless ratio of intracapillary to extracapillary transport resistance from Sharan et al; D*alpha / k*Rc = 0.84
			double sigma = params.getVal(PARAM_VAS_SIGMA);
			// radius of capillary
			double Rc = params.getVal(PARAM_VAS_RC);
			double voxel_volume_cm = std::pow(params.getVal(PARAM_VOXEL_SIZE_CM), 3);
			//vasculature density (in cm/ml)
			double Lv = voxel_volume_cm / (std::pow(Rc, 2) * pi);
			// radius of tissue in Krogh Cylinder 
			double Rt = 1 / std::pow(Lv * pi, 0.5);
			double w = Rc / Rt;
			double lambda = 1 - w * w;
			double Kv = 2 * pi * params.getVal(PARAM_O2_DIFFUSIVITY) * (lambda / (sigma * lambda - (2 + lambda) / 4 + (1 / lambda) * std::log(1 / w)));
			double O2_transport_raw = Kv * Lv * (params.getVal(PARAM_VAS_O2_CONC) - O2);
			double O2_transport = O2_transport_raw > 0 ? O2_transport_raw : 0;
			//std::cout << "Lv: " << Lv << ", Rt: " << Rt << ", w: " << w << ", lambda: " << lambda << ", Kv: " << Kv << ", Current O2 level: " << O2 << ", O2 transport rate: " << O2_transport << std::endl;
			//oxygen transport rate is dependent on tumor compartment.
			pVas->setup_chem_source(pVas->get_source_O2(), CHEM_O2, O2_transport);
			
			//Adding immune cell recruitment source to grid
			int neighborhood_cancer_count = neighborhood_cancer_counts(c_vas.x + neighborhood_radius, c_vas.y + neighborhood_radius, c_vas.z + neighborhood_radius);
			double max_num_cell = neighborhood_voxel_counts(c_vas.x + neighborhood_radius, c_vas.y + neighborhood_radius, c_vas.z + neighborhood_radius);
			double local_cancer_ratio = neighborhood_cancer_count / max_num_cell;	

			//std::cout << "cabo: " << _concentration_cabo << ", PARAM_IC50_VEGFR2: " << params.getVal(PARAM_IC50_VEGFR2) << ", PARAM_IC50_AXL: " << params.getVal(PARAM_IC50_AXL)
			//<< ", PARAM_IC50_MET: " << params.getVal(PARAM_IC50_MET) << ", PARAM_LAMBDA_V_CABO: " << params.getVal(PARAM_LAMBDA_V_CABO) << ", R_cabo: " << _R_cabo
			//<< ", vas_factor: " << (1 - params.getVal(PARAM_LAMBDA_V_CABO) * (_concentration_cabo / (_concentration_cabo + params.getVal(PARAM_IC50_VEGFR2))) * _R_cabo) << std::endl;
			
			
			//double ECM_density = p_agGrid->getECMDensity();
			//double ECM_coeff = 1 - ((ECM_density - params.getVal(PARAM_FIB_ECM_BASELINE)) / (params.getVal(PARAM_FIB_ECM_SATURATION) - params.getVal(PARAM_FIB_ECM_BASELINE)));
			//ECM_coeff = (ECM_coeff > 0 ? ECM_coeff : 0);

			//double max_cancer = _sizeX * _sizeY * _sizeZ;
			//double tumor_scaler = std::sqrt(1e5 * max_cancer / (_concentration_cc * params.getVal(PARAM_AVOGADROS)));

			double ECM_density = get_coord_ECM(c_vas);
			double ECM_sat = (ECM_density - params.getVal(PARAM_FIB_ECM_BASELINE)) / (params.getVal(PARAM_FIB_ECM_SATURATION) - params.getVal(PARAM_FIB_ECM_BASELINE));
			double ECM_coeff = std::exp(-1 * ECM_sat);
			ECM_coeff = (ECM_coeff > 0 ? ECM_coeff : 0);
			double ifng_conc = get_chem(c_vas, CHEM_IFN);
			double H_IFNG = ifng_conc / (ifng_conc + params.getVal(PARAM_TEFF_IFN_EC50));
			
			if (local_cancer_ratio < 1) {
				// vas_scaler is the scaling factor accounting for the number of vasculature in the ABM model
				// more vas cell does not mean more T cell entry
				double vas_scaler = 100 / double(get_stats().getVas());
				//double p_entry =  H_IFNG  * tumor_scaler * vas_scaler;
				double p_entry =  H_IFNG  * tumor_scaler * vas_scaler;
				// std::cout << "ECM_coeff: " << ECM_coeff << ", p_entry: " << p_entry << std::endl;
				if (rng.get_unif_01() < p_entry) {
					_t_source.push_back(c_vas);
				}
			}
			//std::cout << "prob T cell source: " << p_Tcell_source << ", prob MDSC source: " << p_MDSC_source << ", prob MAC source: " << p_Mac_source << std::endl;
		}	
	}

	for_each_grid_coord(true, true, true, [&](Coord3D& c) {//Monocyte recruitment

		//Adding immune cell recruitment source to grid
		int neighborhood_cancer_count = neighborhood_cancer_counts(c.x + neighborhood_radius, c.y + neighborhood_radius, c.z + neighborhood_radius);
		double max_num_cell = neighborhood_voxel_counts(c.x + neighborhood_radius, c.y + neighborhood_radius, c.z + neighborhood_radius);
		double local_cancer_ratio = neighborhood_cancer_count / max_num_cell;	
		
		double ccl2_conc = get_chem(c, CHEM_CCL2);
		double H_CCL2 = ccl2_conc / (ccl2_conc + params.getVal(PARAM_MDSC_EC50_CCL2_REC));
		// double H_CCL2 = (H_CCL2_raw > 0.01 ? H_CCL2_raw : 0);
		// double p_MDSC_source = H_CCL2 * (1 - local_cancer_ratio);			
		if (local_cancer_ratio < 1) {
			double p_MDSC_source = H_CCL2;
			//std::cout << "p_MDSC_source: " << p_MDSC_source << ", H_CCL2: " << H_CCL2 << ", local_cancer_ratio: " << local_cancer_ratio << std::endl;
			if (rng.get_unif_01() < p_MDSC_source) {
				_mdsc_source.push_back(c);
			}
		}
		
		// int neighborhood_cancer_count = neighborhood_cancer_counts(c.x + neighborhood_radius, c.y + neighborhood_radius, c.z + neighborhood_radius);
		// double max_num_cell = neighborhood_voxel_counts(c.x + neighborhood_radius, c.y + neighborhood_radius, c.z + neighborhood_radius);
		// double local_cancer_ratio = neighborhood_cancer_count / max_num_cell;	

		//double p_Mac_source = H_CCL2 * (1 - local_cancer_ratio);
		if (local_cancer_ratio < 1) {
			double p_Mac_source = H_CCL2;
			if (rng.get_unif_01() < p_Mac_source) {
				_mac_source.push_back(c);
			}
		}
		
	});
	std::cout << "t cell vector size: " << _t_source.size() << '\n';
	std::cout << "MDSC cell vector size: " << _mdsc_source.size() << '\n';
	std::cout << "MAC cell vector size: " << _mac_source.size() << '\n';
	return;
	
}

void Tumor::reset_immune_source() {
	_t_source.clear();
	_mdsc_source.clear();
	_mac_source.clear();
}

static bool sortNonTCell(CellAgent * ptrCell) { return ptrCell->getType() != AgentTypeEnum::CELL_TYPE_T; }

/*! simulate generaic tumor compartment for one time slice.
	Steps:
	-# process ABM step (cellular scale events)
		-# recruit cells
		-# for each cell,
			-# progress one step and determine consequences.
			-# process consequence of the step.
		-# go through cell list, remove dead cells, add live cells to ABM snapshot stats
		-# refresh MHC
	-# process molecular scale substeps
	-# shuffle cell list.
*/
void Tumor::timeSlice(unsigned long slice) {

	static double dt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
	const double t = slice * dt;
	//static unsigned int nrMolSlice = params.getVal(PARAM_MOL_STEP_PER_SLICE);
	//static double dtMol = dt / nrMolSlice;
	// reset vasculature entry point based on
	update_vas();
	// update ECM every time step
	update_ECM();

	// reset time step stats
	_stats.resetStats();


	// Reset to zeros
	std::fill(_var_abm_to_qsp.begin(), _var_abm_to_qsp.end(), 0);

	//------------------------   Cellular scale events   -----------------------//
	//BOOST_TEST_MESSAGE("slice:" << slice);
	//cout << "slice:" << slice << endl;

	// Recruitment
	// model comparison: find invasive front to recruit T cells and MDSC

	time_slice_recruitment();

	//cout << "after recruitment:" << endl;

	//std::cout << "RNG check (" << slice << ") rec: " << rng.get_unif_01() << std::endl;

	time_slice_movement(t, dt);
	//cout << "after movement:" << endl;

	//std::cout << "RNG check (" << slice << ") move: " << rng.get_unif_01() << std::endl;
	time_slice_state_change(t, dt);

	//cout << "after state:" << endl;
	//std::cout << "RNG check (" << slice << ") state: " << rng.get_unif_01() << std::endl;

	if (_allow_shift_grid && slice % params.getVal(PARAM_SHIFTGRID_INTERVAL) == 0)
	{
		shift_adjust_center();
	}
	//cout << "after shift:" << endl;
	

	//std::cout << "RNG check (" << slice << ") shift: " << rng.get_unif_01() << std::endl;
	time_slice_final_scan();

	//cout << "after clear dead:" << endl;


	// ------------------------  molecular scale events     ---------------------//
	/*
	for (auto p : _vecAgent)
	{
		if (p->getState() == T_CELL_CYT)
		{
			auto pt = dynamic_cast<TCell*>(p);
			if (pt->get_source_IFNg()->get_index() != _chem.get_voxel_idx(pt->getCoord()))
				if (!pt->get_source_IFNg()->is_to_remove())
					std::cout << "location mismatch (IFNg):" << pt->getCoord() << ", "
						<< _chem.idx_to_coords(pt->get_source_IFNg()->get_index())
						<< std::endl;
		}
	}*/
	time_slice_molecular(t);

	//std::cout << "RNG check (" << slice << ") mol: " << rng.get_unif_01() << std::endl;

	//cout << "after molecular:" << endl;
	// shuffle cell vector
	if (slice % params.getVal(PARAM_SHUFFLE_CELL_VEC_INTERVAL) == 0)
	{
		std::shuffle(_vecAgent.begin(), _vecAgent.end(), rng.getRNG());
	}


	/*
	std::cout << "Num sources: " 
		<< _chem.get_num_source_sink() << std::endl;
	std::cout << "CELL INFO" << std::endl;
	}*/

	// remove all immune source 
	reset_immune_source();
}

/*! Adjust the center of the tumor compartment
*/
void Tumor::shift_adjust_center(void){
	Coord3D c = get_cam_shift();
	// only shift when distance is larger than c_tol
	int c_tol = 0;
	//std::cout << c << "," << c.length() << std::endl;
	//if (c.length()>c_tol)
	if (c.z != 0)
	{
		//shift_grid(c);
		//std::cout << "Shifting: " << c << "," << c.length() << std::endl;
		Coord3D c_zshift(0, 0, c.z);
		shift_grid(c_zshift);
		std::cout << "Shifting: " << c_zshift << std::endl;
	}
	return;
}

/*!
*/
Coord3D Tumor::get_cam_shift(void) {

//! get center of mass
	unsigned int nr_cancer = 0;
	Coord c(0, 0, 0);

	for (auto p : _vecAgent){
		if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
			c = c + p->getCoord();
			nr_cancer += 1;
		}
	}
	//std::cout << "center: " << _center_target << "CoM: "<< c/nr_cancer << "," << c.length() << std::endl;
	if (nr_cancer) {
		return c / nr_cancer - _center_target;
	}
	else {
		return c;// (0,0,0)
	}
}
/*! Shift content of grid by crd
	In current version, grid size change is not allowed.
	When we need to shift center of grid by (x, y, z), 
	we move the content of the grid instead.
	The following contents need to be relocated:
	# Diffusion grid
		# concentration
		# source/sink
	# Agent grid
		# cells
		# structures
			# recruitment sources
	\param [in] crd: shifting distance (x, y, z) (of "camera")
*/
void Tumor::shift_grid(Coord3D& c_shift) {
	// k_step = 1 if k >0 else -1 (direction)
	bool x_pos = c_shift.x >= 0;
	bool y_pos = c_shift.y >= 0;
	bool z_pos = c_shift.z >= 0;

				
	int nr_shift = _sizeX*_sizeY*_sizeZ
		- (_sizeX - abs(c_shift.x))*(_sizeY - abs(c_shift.y))*(_sizeZ - abs(c_shift.z));
	int i_shift = 0;
	//cout << "shift: " << c_shift << ", nr: " << nr_shift << endl;
	CoordVec drop_out(nr_shift, Coord(0, 0, 0));
	/* shift agent grid:
		move camera by (x, y, z): 
		agents: (x0, y0, z0) -> (x0-x, y0-y, z-z0)
		*/
	// first round: save voxels to be moved out of sight to temp grid
	//	and shift remaining voxels to new location 

	for_each_grid_coord(x_pos, y_pos, z_pos, [&](Coord3D& c){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		//BOOST_TEST_MESSAGE(c);
		if (!_agGrid.inGrid(c_new)){
			//add to temp grid
			_agGrid_temp(c_target) = _agGrid(c);
			drop_out[i_shift] = c;
			i_shift++;
		}
		else{
			_agGrid(c_target) = _agGrid(c);
			for (auto pAg : _agGrid(c_target)->get_agents()){
				//BOOST_TEST_MESSAGE(c_target);
				auto pCell = dynamic_cast<Cell_Tumor*>(pAg);
				pCell->setCoord(c_target);
			}
		}
		return;
	});

	//cout << "added to drop_out: " << i_shift << endl;

	// populate new voxels (recycling out-of-site ones)
	/*
	for_each_voxel(x_pos, y_pos, z_pos, [&](Coord3D& c){
		return;
	});*/

	for (auto& c : drop_out){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		if (!_agGrid.inGrid(c_new)){
			// overwrite info and put back to grid
			_agGrid(c_target) = _agGrid_temp(c_target);
			// apply changes to voxels entering the grid (recycled)
			// remove old, 
			for (auto pAg : _agGrid(c_target)->get_agents()){
				auto pCell = dynamic_cast<Cell_Tumor*>(pAg);
				pCell->set_drop_out();
			}
			_agGrid(c_target)->remove_all_agents();
			//populate new
			populate_voxel_random(c_target);
		}
	}

	// remove dropped cells
	// now done together with dead cell

	// shift recruitment sources? 

	/* shift diffusion grid
		# shift values
		# shift point source/sink
		*/
	unsigned int num_sub = _chem.get_num_substrates();
	BioFVMGrid::chem_grid _chem_temp = _chem.get_concentrations();
	for_each_grid_coord(x_pos, y_pos, z_pos, [&](Coord3D& c){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		unsigned int idx = _chem.get_voxel_idx(c_target);
		for (size_t i = 0; i < num_sub; i++)
		{
			if (_agGrid.inGrid(c_new)){
				_chem_temp[idx][i] = _chem(c, i);
			}
			else{
				_chem_temp[idx][i] = 0;
			}
		}
		return;
	});
	_chem.setup_concentrations(_chem_temp);

	// shift source/sink: already taken care of.
	// Sources moving out of grid: set removed with agents
	// Other sources: moved together with agents.
	/* old code
	for (auto pS : _chem.get_sink_source()){
		auto c = _chem.idx_to_coords(pS->get_index());
		auto c_target = c - c_shift;
		_chem.move_source_sink(pS, c_target);
	}*/

	return;
}

unsigned int Tumor::make_unique_id()
{
	return _next_unique_id++;
}


void Tumor::time_slice_recruitment(){


	double T_cell_recruitment_rate = params.getVal(PARAM_TEFF_RECRUIT_K);
	double pRecTeff = get_T_recruitment_prob(_concentration_t_cyt, T_cell_recruitment_rate);

	// CD8
	const auto dummyTcell = _tInitDummy;
	std::cout << "T cell concentration:" << _concentration_t_cyt << "\n"
		<< "recruitment rate: " << params.getVal(PARAM_TEFF_RECRUIT_K) << "\n"
		<< "rec prob (Teff): " << pRecTeff << std::endl;
	for (auto& crd : _t_source){
		if (rng.get_unif_01() < pRecTeff) {
			bool rec = recruitOneCellInMooreNeighborhood(dummyTcell, crd, rng);
			if (rec)
			{
				auto pT = dynamic_cast<TCell*>(_vecAgent.back());
				inc_abm_var_exchange(TUMEX_TEFF_REC);
				_stats.incRecruit(dummyTcell->getType(), dummyTcell->getState());
			}
		}
	}

	// Add Treg recruitment
	double Treg_recruitment_rate = params.getVal(PARAM_TREG_RECRUIT_K);
	double pRecTreg = get_T_recruitment_prob(_concentration_t_reg, Treg_recruitment_rate);

	std::cout << "Treg concentration:" << _concentration_t_reg << "\n"
		<< "recruitment rate: " << params.getVal(PARAM_TREG_RECRUIT_K) << "\n"
		<< "rec prob (Treg): " << pRecTreg << std::endl;
	// Treg recruitment loop
	const auto dummyTreg = _tcd4InitDummy;
	for (auto& crd : _t_source) {
		if (rng.get_unif_01() < pRecTreg) {
			bool rec = recruitOneCellInMooreNeighborhood(dummyTreg, crd, rng);
			if (rec) {
				auto pT = dynamic_cast<TCD4*>(_vecAgent.back());
				pT->setTreg();
				//pT->setup_chem_source(pT->get_source_IFNg(), CHEM_IFN, 0);
				pT->setup_chem_source(pT->get_source_IL_2(), CHEM_IL_2, 0);
				pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, params.getVal(PARAM_TREG_IL_10_RELEASE));
				pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, params.getVal(PARAM_TREG_TGFB_RELEASE));
				//std::cout << "Treg divide interval: " << pT->get_divide_cd_TCD4_exp() << std::endl;
				_stats.incRecruit(dummyTreg->getType(), dummyTreg->getState());
				inc_abm_var_exchange(TUMEX_TREG_REC);
			}
		}
	}
	
	// Add Th recruitment, the default state of TCD4 is Th
	double Th_recruitment_rate = params.getVal(PARAM_TH_RECRUIT_K);
	double pRecTh = get_T_recruitment_prob(_concentration_t_h, Th_recruitment_rate);

	std::cout << "Th concentration:" << _concentration_t_h << "\n"
		<< "recruitment rate: " << params.getVal(PARAM_TH_RECRUIT_K) << "\n"
		<< "rec prob (Th): " << pRecTh << std::endl;
	const auto dummyTh = _tcd4InitDummy;
	// Th recruitment loop
	for (auto& crd : _t_source) {
		if (rng.get_unif_01() < pRecTh) {
			bool rec = recruitOneCellInMooreNeighborhood(dummyTh, crd, rng);
			if (rec) {
				auto pT = dynamic_cast<TCD4*>(_vecAgent.back());
				//pT->setup_chem_source(pT->get_source_IFNg(), CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
				pT->setup_chem_source(pT->get_source_IL_2(), CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
				pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, 0);
				pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, 0);
				//std::cout << "TCD4 divide interval: " << pT->get_divide_cd_TCD4_exp() << std::endl;

				_stats.incRecruit(dummyTh->getType(), dummyTh->getState());
				inc_abm_var_exchange(TUMEX_TH_REC);
			}
		}
	}


	//std::cout << "number of th recruited: " << cd4_recruited << std::endl;

	// MDSC
	/**/
	const auto MDSC_dummy = _MDSCInitDummy;
	//std::cout << "rec prob (MDSC): " << pRec_MDSC << std::endl;	

	const double voxel_volume = std::pow(double(params.getVal(PARAM_VOXEL_SIZE)) / 1e6, 3);

	for (auto& crd : _mdsc_source){	
		double pRec_MDSC = params.getVal(PARAM_MDSC_RECRUIT_K);
		//std::cout << "PARAM_MDSC_RECRUIT_K: " << params.getVal(PARAM_MDSC_RECRUIT_K) << std::endl;


		if (rng.get_unif_01() < pRec_MDSC) {
			bool rec = recruitOneCellInMooreNeighborhood(MDSC_dummy, crd, rng);		
			if (rec)
			{
				_stats.incRecruit(MDSC_dummy->getType(), MDSC_dummy->getState());
			}
		}
	}	

	const auto mac_dummy = _macInitDummy;

	for (auto& crd : _mac_source) {
		double pRec_Mac = params.getVal(PARAM_MAC_RECRUIT_K);

		//std::cout << "PARAM_MAC_RECRUIT_K: " << params.getVal(PARAM_MAC_RECRUIT_K) << std::endl;

		if (rng.get_unif_01() < pRec_Mac) {
			bool rec = recruitOneCellInMooreNeighborhood(mac_dummy, crd, rng);
			if (rec)
			{
				_stats.incRecruit(mac_dummy->getType(), mac_dummy->getState());
				auto pT = dynamic_cast<Mac*>(_vecAgent.back());
				if (pT->getState() == AgentStateEnum::MAC_M1)
				{
					if (rng.get_unif_01() < 0.3) {
						pT->setM2();
					}
				}
			}
		}
	}
}

void Tumor::time_slice_movement(double t, double dt) {
    for (CellVec::size_type i = 0; i < _vecAgent.size(); i++) {
        auto ptP = dynamic_cast<Cell_Tumor*>(_vecAgent[i]);
        
        // For fibroblasts, attempt multiple moves per timestep
        int move_attempts = 1;
        if (auto fib = dynamic_cast<Fib*>(ptP) && t > 0) {
            move_attempts = int(params.getVal(PARAM_FIB_MOVE_STEPS)); 
			//move_attempts = 10;
			//std::cout << "fib move attempts: " << move_attempts << std::endl;
        }
		else if (auto cd4 = dynamic_cast<TCD4*>(ptP) && t > 0) {
			move_attempts = int(params.getVal(PARAM_T_CELL_MOVE_STEPS));
			//move_attempts = 30;
			//std::cout << "cd4 move attempts: " << move_attempts << std::endl;
		}
		else if (auto cd8 = dynamic_cast<TCell*>(ptP) && t > 0) {
			move_attempts = int(params.getVal(PARAM_T_CELL_MOVE_STEPS));
			//move_attempts = 30;
			//std::cout << "cd8 move attempts: " << move_attempts << std::endl;
		}
		else if (auto cancer = dynamic_cast<CancerCell*>(ptP) && t > 0) {
			move_attempts = ptP->getState() == CANCER_STEM ?
				int(params.getVal(PARAM_CANCER_STEM_MOVE_STEPS)) :
				int(params.getVal(PARAM_CANCER_CELL_MOVE_STEPS));
			//move_attempts = 1;
			//std::cout << "cancer move attempts: " << move_attempts << std::endl;
		}
		else if (auto macrophage = dynamic_cast<Mac*>(ptP) && t > 0) {
			move_attempts = int(params.getVal(PARAM_MAC_MOVE_STEPS));
			//move_attempts = 5;
			//std::cout << "macrophage move attempts: " << move_attempts << std::endl;
		}
		else if (auto mdsc = dynamic_cast<MDSC*>(ptP) && t > 0) {
			move_attempts = int(params.getVal(PARAM_MDSC_MOVE_STEPS));
			//move_attempts = 5;
			//std::cout << "mdsc move attempts: " << move_attempts << std::endl;
		}

        for (int attempt = 0; attempt < move_attempts; attempt++) {
            Coord c(0, 0, 0);
            bool move = ptP->agent_movement_step(t, dt, c);
            
            if (move) {
                // Store the initial position
                Coord3D start_pos = ptP->getCoord();
                
                // Remove from initial position
                try {
                    removeAgentFromGrid(start_pos, ptP);
                }
                catch (...) {
                    cout << ptP->toString() << endl;
                    cout << start_pos << "->" << c << endl;
                    std::rethrow_exception(std::current_exception());
                }
                
                // Move to final position
                ptP->setCoord(c);
                addAgentToGrid(c, ptP);

                // Move the tail for fibroblasts
                if (auto fib = dynamic_cast<Fib*>(ptP)) {
                    for (auto tail = fib->getNext(); tail != nullptr; tail = tail->getNext()) {
                        Coord3D start_pos2 = tail->getCoord();
                        removeAgentFromGrid(tail->getCoord(), tail);
                        tail->setCoord(start_pos);
                        addAgentToGrid(start_pos, tail);
                        start_pos = start_pos2;
                    }
                }

                _stats.incMove(ptP->getType(), ptP->getState());
            }
        }
    }
    return;
}

void Tumor::time_slice_state_change(double t, double dt){

	/*Scan neighbor*/
	for (auto const p: _vecAgent)
	{
		auto pCell = dynamic_cast<Cell_Tumor*>(p);
		pCell->agent_state_scan();
	}
	
	//std::cout << "RNG check "<< " scan: " << rng.get_unif_01() << std::endl;
	/*process state change*/
	CellVec::size_type originalCount = _vecAgent.size();
	for (CellVec::size_type i = 0; i < originalCount; i++)
	{
		CellAgent *ptP = _vecAgent[i];
		Coord c(0, 0, 0);

		//std::cout << "processing: " << ptP->getID() << ", " << ptP->getType() << ", " << ptP->getState() << std::endl;
		bool divide = ptP->agent_state_step(t, dt, c);
		//std::cout << "divide: " << divide << std::endl;
		if (divide)
		{
			if (ptP->getType() != AgentTypeEnum::CELL_TYPE_FIB)
			{
				//std::cout << "adding daughter cell" << std::endl;
				auto ptD = ptP->createCellCopy();

				if (ptP->getType() == AgentTypeEnum::CELL_TYPE_CANCER)
				{
					//std::cout << "cancer cell" << std::endl;
					auto ptPCancer = dynamic_cast<CancerCell*>(ptP);
					auto ptDCancer = dynamic_cast<CancerCell*>(ptD);

					ptDCancer->_stemID = ptPCancer->_stemID;
					if (ptP->getState() == AgentStateEnum::CANCER_STEM)
					{
						bool asymmetric = rng.get_unif_01() < params.getVal(PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB);
						if (asymmetric)
						{
							ptDCancer->setProgenitor();
						}
						else {
							ptDCancer->_stemID = ptPCancer->getID();
						}
					}
				}
				if (ptP->getType() == AgentTypeEnum::CELL_TYPE_VAS)
				{
					//std::cout << "cancer cell" << std::endl;
					auto ptPVas = dynamic_cast<Vas*>(ptP);
					auto ptDVas = dynamic_cast<Vas*>(ptD);

					/*
					if (ptPVas->getState() == AgentStateEnum::VAS_STALK && ptPVas->getNext() == nullptr) 
					{
						ptPVas->setNext(ptDVas);
						ptDVas->setPrevious(ptPVas);
						//Need to think about what happened when branching
					}
					*/
					if (ptPVas->getState() == AgentStateEnum::VAS_TIP) {
						//std::cout << "Tip cell parent: " << ptPVas->getCoord() << "Tip cell daughter: " << ptDVas->getCoord() << std::endl;
						//std::cout << "check vas availablility :  success: " << ptPVas->getCoord() << std::endl;
						ptDVas->set_phalanx();
						ptPVas->set_upstream_neighbor(ptDVas);
					}
					else {
						//std::cout << "cell: " << ptDVas->getCoord() << " branched" << std::endl;
						//std::cout << "phalanx cell parent: " << ptPVas->getCoord() << "phalanx cell daughter: " << ptDVas->getCoord() << std::endl;
						ptDVas->set_tip();
						ptDVas->set_upstream_neighbor(ptPVas);
						ptDVas->set_tip_id(make_unique_id());
					}
				
				}

				//add daugther cell
				ptD->setCoord(c);
				//std::cout << "Vas daughter cell: "<< ptD->getCoord() << std::endl;
				addAgentToGrid(c, ptD);
				_vecAgent.push_back(ptD);

				if (ptP->getState() == AgentStateEnum::VAS_TIP) {
					//std::cout << "Tip cell parent: " << ptP->getCoord() << "Tip cell daughter: " << ptD->getCoord() << std::endl;
				}
				if (ptP->getState() == AgentStateEnum::VAS_PHALANX)
				{
					//std::cout << "phalanx cell parent: " << ptP->getCoord() << "phalanx cell daughter: " << ptD->getCoord() << std::endl;
				}

				//auto pD_tumor = dynamic_cast<Cell_Tumor*>(ptD);
				//pD_tumor->move_all_source_sink();

				_stats.incProlif(ptD->getType(), ptD->getState());
				//std::cout << "added" << std::endl;
			}
			if (ptP->getType() == AgentTypeEnum::CELL_TYPE_FIB)
			{
				Fib* ptFib = dynamic_cast<Fib*>(ptP);
				/* two steps :
				   1. connect the fibroblast tail with additional fibroblast agent
				   2. loop from the original tail to the head and change their state
				 */
				if (ptFib->getNext() == nullptr)
				{
					Coord3D c1, c2;
					if (getOpenNeighborForFib(ptFib->getCoord(), c1)) {
						if (getOpenNeighborForFib(c1, c2)) {
							if (c2 != ptFib->getCoord()) {
								// 1. connect the fibroblast tail with additional fibroblast agent
								ptFib->setCAF();
								CellAgent* p1 = createOneInitCell(AgentTypeEnum::CELL_TYPE_FIB, AgentStateEnum::FIB_CAF, c1);
								CellAgent* p2 = createOneInitCell(AgentTypeEnum::CELL_TYPE_FIB, AgentStateEnum::FIB_CAF, c2);
								auto ptF1 = dynamic_cast<Fib*>(p1);
								auto ptF2 = dynamic_cast<Fib*>(p2);

								ptFib->setNext(ptF1);
								ptF1->setPrevious(ptFib);
								ptF1->setNext(ptF2);
								ptF2->setPrevious(ptF1);

								// 2. loop from the original tail to the head and change their state
								for (auto head = ptFib->getPrevious(); head != nullptr; head = head->getPrevious()) {
									head->setCAF();
								}
							}

						}
					}
				}
			}
		}
	}
	return;
}
/*! Final scan during a slice
	Process cell counts and slice statistics 
	Remove dead cells from cell vector
*/
void Tumor::time_slice_final_scan(void){

	CellVec::iterator last = stable_partition(_vecAgent.begin(), _vecAgent.end(), [&](CellAgent* p){
		auto pC = dynamic_cast<Cell_Tumor*>(p);
		return !(pC->isDead() || pC->is_drop_out());
	});// sortLiveCellAddToStats);

	// count live cells (and cancer cells)
	int nrPDL1 = 0;
	for (auto it = _vecAgent.begin(); it < last; it++)
	{
		Cell_Tumor* ptrCell = dynamic_cast<Cell_Tumor*>(*it);
		//statistics
		_stats.incCellState(ptrCell->getType(), ptrCell->getState());
		if (ptrCell->getType() == CELL_TYPE_CANCER)
		{
			inc_abm_var_exchange(TUMEX_CC);
		}
		if (ptrCell->is_PDL1_pos())
		{
			nrPDL1++;
		}
	}
	//update cancer cell counts
	_cancer_counts = get_stats().getCancerCell();
	_stats.set_stats_misc(STATS_MISC_PDL1_POS, nrPDL1);
	int nr_death = 0;
	int nr_dropout = 0;
	// delete removed cells
	for (CellVec::iterator it = last; it < _vecAgent.end(); it++)
	{
		auto p = dynamic_cast<Cell_Tumor*>(*it);
		if (p->isDead())
		{

			// Cancer cell dies -> antigen presentation
			if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
				inc_abm_var_exchange(TUMEX_CC_DEATH);
				//Release antigen after death
			}
			_stats.incDeath(p->getType(), p->getState());
			if (!p->is_drop_out())
			{
				// drop out cells are already removed from voxel
				removeAgentFromGrid(p->getCoord(), p);
			}
			nr_death++;
		}
		else if (p->is_drop_out() && !p->isDead()){
			// dead and dropout are only counted as dead

			removeAgentFromGrid(p->getCoord(), p);
			_stats.incDropOut(p->getType(), p->getState());
			nr_dropout++;
		}
		else{
			std::cout << "Error: cells neither dead nor drop out are removed" << std::endl;
		}
		// Remove from grid when step is finished so they can still function during the state_change step
		delete p;
	}
	_vecAgent.erase(last, _vecAgent.end());
	//std::cout << "Death: " << nr_death << ", dropout: " << nr_dropout << std::endl;
	return;
}

void Tumor::time_slice_molecular(double t){

	static double dt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
	_chem.remove_dead_source_sink();
	if (params.getVal(PARAM_DIFFUSION_ON) && !_chem.grid_skippable())
	{
		for (auto const pS : _chem.get_sink_source()){
			auto c = _chem.idx_to_coords(pS->get_index());
		}
		static double dt_mol = dt / params.getVal(PARAM_MOLECULAR_STEP_PER_SLICE);
		double t_mol = t;
		double t_mol_end = t + dt;
		while (t_mol < t_mol_end)
		{
			// diffusion time step, including decay/point source/sink
			_chem.timestep(dt_mol);
			t_mol += dt_mol;
		}
	}
}




//! default header for extra remark column when writing cell grid to file
std::string Tumor::getExtraRemarkHeader() const {
	std::stringstream ss;
	ss << "extra";
	return ss.str();
}

/*! Print grid snapshot to a snapshot file.
	Info includes: Element type, IL2 concentration, and MHC concentration.
*/
std::string Tumor::printGridToFile() const {
//void Tumor::printGridToFile(unsigned long slice) const {
	std::stringstream ss;
	// header
	ss << "x,y,z,"<<_chem.get_substrate_names() << std::endl;
	// content
	ss << _chem;
	return ss.str();
}

/*! print gradient of cytokines along selected direction
* x=0; y=1; z=2 (axis)
*/
std::string Tumor::printGradientToFile(int axis) {
	//void Tumor::printGridToFile(unsigned long slice) const {
	
	std::stringstream ss;
	// header
	ss << "x,y,z," << _chem.get_substrate_names() << std::endl;
	// content
	for (size_t i = 0; i < _chem.get_num_voxel(); i++) {
		Coord3D c = _chem.idx_to_coords(i);
		std::vector< std::vector<double> > gradient = _chem.get_gradient(i);
		ss << c.x << ',' << c.y << ',' << c.z << ",";
		auto delim = "";
		for (size_t j = 0; j < _chem.get_num_substrates(); j++) {
			ss << delim << gradient[j][axis];
			delim = ",";
		}
		ss << std::endl;
	}
	
	return ss.str();
}

// Ask this question about const (current cannot add const to the function
std::string Tumor::printStromaToFile() {
	//void Tumor::printGridToFile(unsigned long slice) const {
	std::stringstream ss;
	// header
	ss << "x,y,z,ECM_density,Fib_field" << std::endl;
	// content
	
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {

		TumorGridVoxel* p_agGrid = dynamic_cast<TumorGridVoxel*>(_agGrid(c));
		double ECM_density = p_agGrid->getECMDensity();
		double Fib_field = p_agGrid->getFibField();
		//std::cout << c.x << "," << c.y << "," << c.z << "," << density << std::endl;
		ss << c.x << "," << c.y << "," << c.z << "," << ECM_density <<"," << Fib_field << std::endl;
		return;
	});
	return ss.str();
	
}

/*! Print state of ODE system of all T cells to ODE stats file,
	includeing a time stamp and cell ID.
	\param [in] slice: time slice
*/
void Tumor::printCellOdeToFile(unsigned long slice) const {
}

/*! print grid snapshot to screen
	\param [in] slice: time slice
*/
void Tumor::printGridToScreen(unsigned long slice)const {

	// header
	cout << "Time slice: "<< slice << endl;
	cout << "Grid size: " << getGridSize() << ", "<< _sizeX << " by "
		<< _sizeY << " by " << _sizeY <<  endl;
	//cout << getGridContent() << endl;

	cout << "x, y, z, nr_ag," << endl;
	for (int k = 0; k < _sizeZ; k++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int i = 0; i < _sizeX; i++)
			{
				Coord c(i, j, k);
				cout << c.x << "," << c.y << "," << c.z << ",";
				AgentGridVoxel* voxel = _agGrid.get(i, j, k);
				TumorGridVoxel * v2 = dynamic_cast<TumorGridVoxel*>(voxel);
				cout << voxel->getNumAgents();
				//cout << v2->getDistToOrigin();
				//cout << v2->getLoc();
				cout << endl;
			}
		}
	}
}

double Tumor::get_chem(const Coord3D&c, chem_ID i)const{
	if (!_agGrid.inGrid(c))
	{
		return 0;
	}
	else 
	{ 
		return _chem(c, i); 
	}
}

std::vector< std::vector<double> > Tumor::get_gradient_voxel(const Coord3D& c) {
	if (!_agGrid.inGrid(c))
	{
		return { {0} };
	}
	else
	{
		size_t idx = _chem.get_voxel_idx(c);
		std::cout << "coordinate index: " << idx << std::endl;
		return _chem.get_gradient(idx);
	}
}
//-----------------------------------------  Protected --------------------------------------//

//-----------------------------------------  Private  --------------------------------------//

/*! Setup compartment environment. This includes: sources of T cells and MDSCs.
	create vasculature by mapping graph file to the grid.
	all locations mapped onto the graph are designated as T cell or MDSC sources.
*/
void Tumor::initEnvironment() {

	int nrSource = 0;

}

/*! Setup initial cells. Called when simulation starts.
	Create initial cells and set their coordinates;
	put initial cells on the grid;
	put initial cells into cell vector.
	After all intial cells are generated, go through the vector and record initial stats
	\param [in] initialCellFileName: initial arrangement of cells

	read initial cell configuration from a separate file
	structure of initial condition file:
	<initialCondition>
	  <cell>
		<x></x>
		<y></y>
		<z></z>
		<celltype></celltype>
		<cellstate></cellstate>
	  </cell>
	  <cell>
		...
	  </cell>
	</initialCondition>

	<cellCluster/> elements will have an additional subelement <count/>,
	indicating the number of cells in that cluster.
*/
void Tumor::initCell(string initialCellFileName){

	namespace pt = boost::property_tree;
	const std::string cellPath = "initialCondition";
	const std::string cellTag = "cell";
	const std::string cellClusterTag = "cellCluster";
	using std::pair;
	pt::ptree tree;

	// in case no initial condition file provided
	if (initialCellFileName.empty())
	{
		init_cell_fill_grid();

		//init_cell_single_center();
	}
	else {
		try {
			pt::read_xml(initialCellFileName, tree, pt::xml_parser::trim_whitespace);
			// get nodes
			BOOST_FOREACH(pt::ptree::value_type const& cell, tree.get_child(cellPath)) {
				if (cell.first == cellTag) {
					icProperty ic = cell.second;
					unsigned int e = ic.get<unsigned int>("celltype");
					unsigned int s = ic.get<unsigned int>("cellstate", 0);
					unsigned int x = ic.get<unsigned int>("x");
					unsigned int y = ic.get<unsigned int>("y");
					unsigned int z = ic.get<unsigned int>("z");
					AgentType type = static_cast<AgentType>(e);
					AgentState state = static_cast<AgentState>(s);
					createOneInitCell(type, state, Coord3D(x,y,z));
				}
				else if (cell.first == cellClusterTag)
				{
					icProperty ic = cell.second;
					createClusterInitCell(ic);
				}
			}
		}
		catch (std::exception & e) {
			std::cerr << "Error creating initial cells" << std::endl;
			//std::cerr << e.what() << std::endl;
			throw e;
		}
	}

	// check and update stats
	for (auto ptrCell : _vecAgent)
	{
		_stats.incCellState(ptrCell->getType(), ptrCell->getState());
		if (ptrCell->isDead())
		{
			cerr << "dead cell initiated" << endl;
			exit(1);
		}
	}
}

/*! 
*/
void Tumor::init_cell_single_center(void){
	auto crd = Coord3D(_sizeX, _sizeY, _sizeZ) / 2;
	AgentType type = AgentTypeEnum::CELL_TYPE_CANCER;
	AgentState state = AgentStateEnum::CANCER_STEM;
	createOneInitCell(type, state, crd);
	//_stats.incRecruit(type, state);
	return;
}

bool Tumor::getOpenNeighborForFib(Coord3D start, Coord3D& result) {
	/* Two fibroblast agents must be in the vonNeumann neighborhood*/
	const auto shape = _fibInitDummy->getCellShape();
	int idx;
	if (getOneOpenVoxel(shape->getFibDestinationVoxels(),
		shape->getFibDirectionAnchor(), start, AgentTypeEnum::CELL_TYPE_FIB, idx, rng))
	{
		result = shape->getFibDirectionAnchor()[idx] + start;
		return true;
	}
	return false;
}

/*! initial cell: fill grid with random population
*/
void Tumor::init_cell_fill_grid(void){
	int nr_voxel = _agGrid.getSize();
	//std::cout << "start population" << std::endl;
	Coord c_sum(0, 0, 0);
	unsigned int n = 0;
	/*
	// For testing
	unsigned int k = 100;
	std::vector<int> voxel_id(nr_voxel);
	std::iota(std::begin(voxel_id), std::end(voxel_id), 0); 
	//std::shuffle(voxel_id.begin(), voxel_id.begin()+k, rng.getRNG());
	rng.shuffle_first_k<int>(voxel_id, k);
	for (size_t i = 0; i < k; i++)
	{
		Coord3D crd = _agGrid.get_coord(voxel_id[i]);
		//std::cout << "Voxel_id: " << voxel_id[i] << ", " << crd << endl;
		createOneInitCell(AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_CYT, crd);
	}*/ 

	/* For initial condition, cancer cells are generated first as they occur first.
	* Then, fibroblast is generated at the second loop (for_each_grid_coord)
	*/
	for_each_grid_coord(true, true, true, [&](Coord3D& c){
		//initialize Cancer cells
		if (populate_voxel_random(c))
		{
			c_sum = c_sum + c;
			n += 1;
		}
		return;
	});
	
	//std::cout << "_concentration_cc: " << _concentration_cc * params.getVal(PARAM_AVOGADROS) << ", _cmax: " << _cmax << ", _concentration_caf: " << _concentration_caf * params.getVal(PARAM_AVOGADROS)<< std::endl;
	double p_fib_gen = _concentration_caf / (_concentration_cc + _concentration_caf + 1 / params.getVal(PARAM_AVOGADROS));
	if (p_fib_gen > 0.1) {
		p_fib_gen = 0.1;
	}

	if (p_fib_gen < 0.001) {
		p_fib_gen = 0.001;
	}

	std::cout << "p_fib_gen: " << p_fib_gen << std::endl;
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		// Set initial ECM density
		auto p_agGrid = dynamic_cast<TumorGridVoxel*>(_agGrid(c));
		p_agGrid->setECMDensity(params.getVal(PARAM_FIB_ECM_SATURATION));

		Coord3D c1, c2, c3;
		if (rng.get_unif_01() < p_fib_gen) {
			if (getOpenNeighborForFib(c, c1)) {
				if (getOpenNeighborForFib(c1, c2)) {
					if (getOpenNeighborForFib(c2, c3)) {
						if (c3 != c1) {
							//std::cout << c1 << ", " << c2 << ", " << c3 << " " << std::endl;
							CellAgent* p1 = createOneInitCell(AgentTypeEnum::CELL_TYPE_FIB, AgentStateEnum::FIB_NORMAL, c1);
							CellAgent* p2 = createOneInitCell(AgentTypeEnum::CELL_TYPE_FIB, AgentStateEnum::FIB_NORMAL, c2);
							CellAgent* p3 = createOneInitCell(AgentTypeEnum::CELL_TYPE_FIB, AgentStateEnum::FIB_NORMAL, c3);

							Fib* f1 = dynamic_cast<Fib*>(p1);
							Fib* f2 = dynamic_cast<Fib*>(p2);
							Fib* f3 = dynamic_cast<Fib*>(p3);

							f1->setNext(f2);
							f2->setPrevious(f1);
							f2->setNext(f3);
							f3->setPrevious(f2);

						}
					}
				}
			}
		}
	});

	double p_tcell_gen = 0 * (_concentration_t_eff_tum / (_concentration_cc + _concentration_t_eff_tum + 1 / params.getVal(PARAM_AVOGADROS)) );
	std::cout << "p_tcell_gen: " << p_tcell_gen << std::endl;
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		if (rng.get_unif_01() < p_tcell_gen) {
			CellAgent* ptrCell = createOneInitCell(AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_EFF, c);
			if (ptrCell) {
				auto pT = dynamic_cast<TCD4*>(ptrCell);
			}
		}
	});

	double p_cd4_gen =  0.03 * (_concentration_t_h_tum / (_concentration_cc + _concentration_t_h_tum + 1 / params.getVal(PARAM_AVOGADROS)) );
	std::cout << "p_cd4_gen: " << p_cd4_gen << std::endl;
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		if (rng.get_unif_01() < p_cd4_gen) {
			CellAgent* ptrCell = createOneInitCell(AgentTypeEnum::CELL_TYPE_TCD4, AgentStateEnum::TCD4_Th, c);
			if (ptrCell) {
				auto pT = dynamic_cast<TCD4*>(ptrCell);
				//pT->setup_chem_source(pT->get_source_IFNg(), CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
				pT->setup_chem_source(pT->get_source_IL_2(), CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
				pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, 0);
				pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, 0);
			}
		}
	});

	double p_treg_gen = 0 * (_concentration_t_reg_tum / (_concentration_cc + _concentration_t_reg_tum + 1 / params.getVal(PARAM_AVOGADROS)) );
	std::cout << "p_treg_gen: " << p_treg_gen << std::endl;
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		if (rng.get_unif_01() < p_treg_gen) {
			CellAgent* ptrCell = createOneInitCell(AgentTypeEnum::CELL_TYPE_TCD4, AgentStateEnum::TCD4_TREG, c);
			if (ptrCell) {
				auto pT = dynamic_cast<TCD4*>(ptrCell);
				pT->setup_chem_source(pT->get_source_IL_2(), CHEM_IL_2, 0);
				pT->setup_chem_source(pT->get_source_IL_10(), CHEM_IL_10, params.getVal(PARAM_TREG_IL_10_RELEASE));
				pT->setup_chem_source(pT->get_source_TGFB(), CHEM_TGFB, params.getVal(PARAM_TREG_TGFB_RELEASE));
			}
		}
	});

	double p_mdsc_gen = (_concentration_mdsc / (_concentration_cc + _concentration_mdsc + 1 / params.getVal(PARAM_AVOGADROS)) );
	std::cout << "p_mdsc_gen: " << p_mdsc_gen << std::endl;
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		if (rng.get_unif_01() < p_mdsc_gen) {
			CellAgent* ptrCell = createOneInitCell(AgentTypeEnum::CELL_TYPE_MDSC, AgentStateEnum::DEFAULT_STATE, c);
		}
	});

	double p_mac_gen = (_concentration_m2 / (_concentration_cc + _concentration_m2 + 1 / params.getVal(PARAM_AVOGADROS)) );
	std::cout << "p_mac_gen: " << p_mac_gen << std::endl;
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		if (rng.get_unif_01() < p_mac_gen) {
			CellAgent* ptrCell = createOneInitCell(AgentTypeEnum::CELL_TYPE_MAC, AgentStateEnum::MAC_M2, c);
		}
	});

	for (const Coord3D& c : _init_vas) {
		//std::cout << "vas location: (" << c.x << ", " << c.y << ", " << c.z << ") " << std::endl;
		CellAgent* ptrCell = createOneInitCell(AgentTypeEnum::CELL_TYPE_VAS, AgentStateEnum::VAS_PHALANX, c);
		if (ptrCell) {
			//std::cout << "initial vas location: (" << c.x << ", " << c.y << ", " << c.z << ") " << std::endl;
			Vas* pVas = dynamic_cast<Vas*>(ptrCell);
			pVas->set_moveDirection();
			double p_branch = params.getVal(PARAM_VAS_BRANCH_PROB)/5;
			if (rng.get_unif_01() < p_branch) {
				//std::cout << "branch coord: " << pVas->getCoord() << std::endl;
				pVas->set_branch();
			}

			c_sum = c_sum + c;
			n += 1;
		}
	}
	
	std::cout << "finished vas generation" << std::endl;
	_center_target = n ? c_sum / n : c_sum;
	//std::cout << _center_target << std::endl;
	return;
}

Cell_Tumor* Tumor::populate_voxel_random(const Coord3D& crd) {
	AgentType type = AgentTypeEnum::AGENT_DUMMY;
	AgentState state = AgentStateEnum::DEFAULT_STATE;
	//bool create_cancer_cell = false;
	Cell_Tumor* pTumorCell = NULL;
	int div = 0;
	//Adding cells based on initial condition
	if (_voxel_ic.get_type_state(crd, rng, type, state, div)) {
		CellAgent* ptrCell = createOneInitCell(type, state, crd);
		if (ptrCell && type == AgentTypeEnum::CELL_TYPE_CANCER)
		{
			//create_cancer_cell = true;
			CancerCell* pCancerCell = dynamic_cast<CancerCell*>(ptrCell);
			pTumorCell = pCancerCell;
			if (state == AgentStateEnum::CANCER_PROGENITOR)
			{
				pCancerCell->setDivCounter(div);
				pCancerCell->randomize_div_cd(int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) + .5));
			}
			else if (state == AgentStateEnum::CANCER_STEM)
			{
				pCancerCell->randomize_div_cd(int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) + .5));
				
			}
		}
	}

	return pTumorCell;
}

/*! recruit a cluster of cells to the lattice
	\param [in] ic: initial condition

	-# create an empty queue; count = 0
	-# push indicated location to the queue
	-# recruit first cell to lattice, count++
	-# iterate until count matched or queue is empty
	  -# deque, get coordinate
	  -# iterate until no space found
		-# push found location to queue
		-# recruit, count++
*/
void Tumor::createClusterInitCell(icProperty &ic) {

	unsigned int nrCellToCreate = ic.get<unsigned int>("count");
	//ElementType e; not working, cannot directly map string in property_tree to enum
	unsigned int e = ic.get<unsigned int>("celltype");
	unsigned int s = ic.get<unsigned int>("cellstate", 0);
	unsigned int x = ic.get<unsigned int>("x");
	unsigned int y = ic.get<unsigned int>("y");
	unsigned int z = ic.get<unsigned int>("z");


	AgentType type = static_cast<AgentType>(e);
	AgentState state = static_cast<AgentState>(s);
	const ShapeBase *shape;
	switch (type)
	{
	case AgentTypeEnum::AGENT_DUMMY:
		break;
	case AgentTypeEnum::CELL_TYPE_CANCER:
		shape = _cInitDummy->getCellShape();
		break;
	case AgentTypeEnum::CELL_TYPE_T:
		shape = _tInitDummy->getCellShape();
		break;	
	case AgentTypeEnum::CELL_TYPE_TCD4:
		shape = _tcd4InitDummy->getCellShape();
		break;
	case AgentTypeEnum::CELL_TYPE_MDSC:
		shape = _MDSCInitDummy->getCellShape();
		break;	
	case AgentTypeEnum::CELL_TYPE_FIB:
		shape = _fibInitDummy->getCellShape();
		break;
	case AgentTypeEnum::CELL_TYPE_VAS:
		shape = _VasInitDummy->getCellShape();
		break;
	default:
		throw std::invalid_argument("unknown cell type in initial cells");
	}

	unsigned int count = 0;
	std::queue<Coord> nextCenter;
	auto crd = Coord(x, y, z);
	nextCenter.push(crd);
	createOneInitCell(type, state, crd);
	count++;

	while (count < nrCellToCreate || nextCenter.empty())
	{
		Coord c = nextCenter.front();
		nextCenter.pop();
		bool spaceFound = true;
		Coord cNew;
		while (count < nrCellToCreate)
		{
			// create one cell
			int idx;
			//spaceFound = getOneOpenDestinationByScan(shape->getProlifNewOccupy(),
			//	shape->getProlifSearchSeq(), shape->getProlifOccupyMap(),
			//	shape->getProlifRelocateMap(), c, ElementType(e), idx);
			spaceFound = getOneOpenVoxel(shape->getProlifDestinationVoxels(),
				shape->getProlifDestinationAnchor(), c, type, idx, rng);
			if (spaceFound)
			{
				cNew = shape->getProlifDestinationAnchor()[idx] + c;
				createOneInitCell(type, state, cNew);
				count++;
				nextCenter.push(cNew);
			}
			else {
				break;
			}
		}
	}
}

/*! recruit one cell to grid
	\param [in] e: cell type to create
	\param [in] state: cell state
	\param [in] crd: 3D coordinate, const
*/
CellAgent* Tumor::createOneInitCell(AgentType type, AgentState state, const Coord3D& crd) {

	CellAgent * ptrCell = NULL;

	//std::cout << "type: " << type << "; state: " << state << std::endl;
	// only add cell if the target voxel can take it
	if (_agGrid(crd)->isOpenToType(type))
	{
		switch (type)
		{
		case AgentTypeEnum::AGENT_DUMMY:
			throw std::invalid_argument("dummy type in initial cells");
			break;
		case AgentTypeEnum::CELL_TYPE_CANCER:
			ptrCell = _cInitDummy->createCellCopy();
			//std::cout << "create cancer cell" << std::endl;
			break;
		case AgentTypeEnum::CELL_TYPE_T:
			ptrCell = _tInitDummy->createCellCopy();
			break;	
		case AgentTypeEnum::CELL_TYPE_TCD4:
			ptrCell = _tcd4InitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_MAC:
			ptrCell = _macInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_FIB:
			ptrCell = _fibInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_MDSC:
			ptrCell = _MDSCInitDummy->createCellCopy();
			break;	
		case AgentTypeEnum::CELL_TYPE_VAS:
			ptrCell = _VasInitDummy->createCellCopy();
			break;
		default:
			throw std::invalid_argument("unknown cell type in initial cells");
		}

		ptrCell->setCoord(crd);
		ptrCell->setAgentState(state);
		addAgentToGrid(crd, ptrCell);
		// reset state
		// now all starts from base state; otherwise report error
		/*
		if (state != AgentStateEnum::DEFAULT_STATE)
		{
		throw std::invalid_argument("invalide initial cell state");
		}
		*/
		// randomize life
		switch (type)
		{
		case AgentTypeEnum::AGENT_DUMMY:
			break;
		case AgentTypeEnum::CELL_TYPE_CANCER:{
			auto pCancerCell = dynamic_cast<CancerCell*>(ptrCell);
			
			//Then, if they need to be set to progenitor cells, they are set by following code
			if (state == AgentStateEnum::CANCER_PROGENITOR)
			{
				pCancerCell->setProgenitor();
		
			}
			else if (state == AgentStateEnum::CANCER_SENESCENT)
			{
				pCancerCell->setSenescent();
			}
			break;
		}
		case AgentTypeEnum::CELL_TYPE_T:{
			auto pTCell = dynamic_cast<TCell*>(ptrCell);
			int life = (int)(rng.get_unif_01() * pTCell->getCellLife() + 0.5);
			pTCell->setCellLife(life);
			break;
		}
		default:
			break;
		}
		// add to cell vector
		_vecAgent.push_back(ptrCell);
		_stats.incDropIn(type, state);
	}

	return ptrCell;
}

/*! return variables needed for QSP module 
	already updated during simulation.
*/
const std::vector<double>& Tumor::get_var_exchange(void){
	return _var_abm_to_qsp;
}
/*! update ABM module with variables from QSP 
*/
void Tumor::update_abm_with_qsp(const std::vector<double>& qsp_var){
	_concentration_cc = qsp_var[LymphCentral::QSPEX_TUM_C] / params.getVal(PARAM_AVOGADROS);
	_concentration_t_cyt = qsp_var[LymphCentral::QSPEX_CENT_TEFF];
	_concentration_t_reg = qsp_var[LymphCentral::QSPEX_CENT_TREG];
	_concentration_t_h = qsp_var[LymphCentral::QSPEX_CENT_TH];
	_concentration_t_eff_tum = qsp_var[LymphCentral::QSPEX_TUM_TEFF];
	_concentration_t_reg_tum = qsp_var[LymphCentral::QSPEX_TUM_TREG];
	_concentration_t_h_tum = qsp_var[LymphCentral::QSPEX_TUM_TH];
	_concentration_m1 = qsp_var[LymphCentral::QSPEX_TUM_M1];
	_concentration_m2 = qsp_var[LymphCentral::QSPEX_TUM_M2];
	_concentration_mdsc = qsp_var[LymphCentral::QSPEX_TUM_MDSC];
	_concentration_nivo = qsp_var[LymphCentral::QSPEX_TUM_NIVO];
	_concentration_ipi = qsp_var[LymphCentral::QSPEX_TUM_IPI];
	_concentration_cabo = qsp_var[LymphCentral::QSPEX_TUM_CABO];
	_concentration_cx = qsp_var[LymphCentral::QSPEX_TUM_CX];
	_concentration_t_exh = qsp_var[LymphCentral::QSPEX_TUM_TEXH];
	_tumor_volume = qsp_var[LymphCentral::QSPEX_TUM_VOL];
	_cmax = qsp_var[LymphCentral::QSPEX_TUM_CMAX];
	_concentration_caf = qsp_var[LymphCentral::QSPEX_TUM_CAF];
	_f_tum_cap = _concentration_cc / _cmax;
	return;
}

/*! T recruitment probability
*/
double Tumor::get_T_recruitment_prob(double c, double k_rec) const{
	/*
	p = k (1/mol) * Cent.T (mol)
	*/	
	double num_rec = c * k_rec;
	//std::cout << "_tumor_volume: " << _tumor_volume << ", scaling_factor: " << scaling_factor << std::endl;
	double p = (num_rec < 1 ? num_rec : 1);

	return p;
}

/*! MDSC recruitment probability
*/
// calculate recruitment probability for both macrophage and MDSCs (Myeloid cells)
double Tumor::get_Myeloid_recruitment_prob(double c, double k_rec) const{
	//p = k (m^3/mol) * Tum.Vol (m^3) * (Tum.MDSCmax * Tum.Vol - Tum.MDSC) (mol)

	//double num_rec = c * k_rec / _tumor_volume;
	double num_rec = c * k_rec;
	//std::cout << "Myeloid Recruitment rate: " << c << ", volume: " << k_rec << ", prob: " << num_rec << std::endl;
	double p;
	if (num_rec > 0){
		p = (num_rec < 1 ? num_rec : 1);
	}
	else{
		p = 0;
	}

	return p;
	
}

void Tumor::update_R_cabo() {
	_R_cabo = _concentration_cabo / (_concentration_cabo + params.getVal(PARAM_IC50_AXL));
	std::cout << "_R_cabo: " << _R_cabo << ", PARAM_IC50_AXL: " << params.getVal(PARAM_IC50_AXL) << std::endl;
}


double Tumor::get_coord_ECM(const Coord3D& c) const {
	TumorGridVoxel* p_agGrid = dynamic_cast<TumorGridVoxel*>(_agGrid(c));
	return p_agGrid->getECMDensity();
}

/*
struct DistancesToNearestQueueEntry {
	Coord3D origin;
	Coord3D coord;

	double distance()const {
		Coord3D diff = coord - origin;
		return std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
	}
};

 A simple KD tree method to find cloest agent given the each coordinate.	
bool operator<(DistancesToNearestQueueEntry const& a, DistancesToNearestQueueEntry const& b) {
	// > instead of < because lower distances are higher priority
	return a.distance() > b.distance();
}

Grid3D<double> Tumor::distances_to_nearest(std::function<bool(Coord3D const&)> callback) {
	std::priority_queue<DistancesToNearestQueueEntry> frontier;
	for_each_grid_coord(true, true, true, [&](Coord3D& c) {
		if ((callback)(c)) {
			frontier.push(DistancesToNearestQueueEntry{ c, c });
		}
		});
	Grid3D<double> distance_to_nearest(_sizeX, _sizeY, _sizeZ, std::numeric_limits<double>::infinity());
	while (!frontier.empty()) {
		DistancesToNearestQueueEntry current = frontier.top();
		frontier.pop();
		if (distance_to_nearest(current.coord) == std::numeric_limits<double>::infinity()) {
			distance_to_nearest(current.coord) = current.distance();

			for (int i = -1; i < 2; i++)
			{
				for (int j = -1; j < 2; j++)
				{
					for (int k = -1; k < 2; k++)
					{
						if (i || j || k) {
							Coord3D neighbor = current.coord + Coord3D(i, j, k);
							if (_agGrid.inGrid(neighbor)) {
								frontier.push(DistancesToNearestQueueEntry{ current.origin, neighbor });
							}
						}
					}
				}
			}
		}
	}

	return distance_to_nearest;
}
*/

};
};
