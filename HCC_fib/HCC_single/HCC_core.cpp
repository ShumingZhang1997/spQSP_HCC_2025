#include "HCC_core.h"

//#include "HCC/SP_QSP_HCC/core/Param.h"
#include "../SP_QSP_HCC/core/GlobalUtilities.h"

#include "SP_QSP_shared/ABM_Base/Grid3D.h"
#include "SP_QSP_shared/ABM_Base/Graph3D.h"

#include "InitialCondition.h"
#include <cmath>         // exp
#include <algorithm>    // std::max
#include <numeric>      // accumulate

extern FileOutputHub output_hub;

extern RNG rng;

//extern SP_QSP_IO::Param params;

namespace SP_QSP_IO {
	namespace SP_QSP_HCC {
		extern Param params;
	}
};
static auto& params = SP_QSP_IO::SP_QSP_HCC::params;

extern InitialCondition ic;
extern std::string initialCellFileName_core;
extern std::string initialCellFileName_margin;

typedef SP_QSP_IO::Coord Coord;

HCC_Core::HCC_Core()
	: _tumor(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y),
		params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z))
	, _lymph()
{


}


HCC_Core::~HCC_Core()
{
}

/*! Setup QSP module.
*/
void HCC_Core::setup_qsp(CancerVCT::Param& p) {
	_lymph.setup_param(p);
	params.update_from_qsp();
}
/*! initialize compartments: randomly populate voxels
	This function is called only when creating the model for the first time,
	and not reading from saved state. Objects created in this function should
	be already in the serialization and can be loaded directly.
*/

void HCC_Core::initializeSimulation(void) {

	// before simulation starts, update tumor with QSP variables
	auto qsp_var = _lymph.get_var_exchange();
	_tumor.update_abm_with_qsp(qsp_var);
	//0.44 is just an arbitary scaler
	//double radius = 0.44 * params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X);
	double radius = 0.44 * params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X);
	std::cout << "presimulation radius: " << radius << std::endl;
	std::cout << "x size: " << params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X) << ", y size: " << params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y) 
			  << ", z size: " << params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z) << std::endl;
	// size = midpoint - factor * size * diameter


	// rule to randomly populate voxel during initlaization or grid shifting
	_tumor.set_allow_shift(ic.getVal(IC_GRID_SHIFT));
	_tumor._voxel_ic.setup(ic.getVal(IC_STATIONARY),
		double(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X) / 2),
		double(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y) / 2),
		double(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z) / 2),
		radius);
	/*
	_tumor._voxel_ic.setup(ic.getVal(IC_DENSITY_CSC),
		int(x_lim),
		int(y_lim),
		int(z_lim),
		int(x_begin),
		int(y_begin),
		int(z_begin));
		*/
		/*
		ic.getVal(IC_X_SIZE),
			ic.getVal(IC_Y_SIZE),
			ic.getVal(IC_Z_SIZE),
			ic.getVal(IC_X_MIN),
			ic.getVal(IC_Y_MIN),
			ic.getVal(IC_Z_MIN)
			*/
			//t cell sources
	{

		using SP_QSP_IO::Graph3D;
		using SP_QSP_IO::Grid3D;

		std::vector<Coord> c_tumor;
		int nr_source_tumor = 0;
		// use vasculature specified in xml file
		if (ic.getVal(IC_USE_FILE_VAS))
		{
			auto vas = Grid3D<int>(params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X),
				params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y),
				params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z), 0);
			Graph3D g;

			auto filename = "../resource/vas.xml";
			g.read_graph(filename);
			g.rasterize_graph(vas, params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_VOXEL_SIZE));

			for (int i = 0; i < vas.getSize(); i++)
			{
				if (vas[i])
				{
					c_tumor.push_back(vas.get_coord(i));
				}
			}

		}
		// randomly choose any location
		else {
			int x_size = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X);
			int y_size = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Y);
			int z_size = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_Z);
			double center_x = x_size / 2.0;
			double center_y = y_size / 2.0;
			double center_z = z_size / 2.0;
			double radius = 0.5 * x_size;
			
			int num_segments = 1;
			for (int seg = 0; seg < num_segments; seg++) {
				Coord current;
				double dx, dy, dz;
				int target_x, target_z;
				// Use z_size - 1 and x_size - 1 to stay within bounds
				double z_coeff = 0.90 + (rng.get_unif_01() * 0.05); 
				double x_coeff = 0.90 + (rng.get_unif_01() * 0.05); 
				
				switch (seg % 4) {
					case 0: // x=0 to x=xmax
						current = Coord(0, center_y, (z_coeff * (z_size - 1)));
						target_x = x_size - 1;
						target_z = current.z;
						dx = 1;
						break;
						
					case 1: // z=0 to z=zmax
						current = Coord((x_coeff * (x_size - 1)), center_y, 0);
						target_x = current.x;
						target_z = z_size - 1;
						dz = 1;
						break;
						
					case 2: // x=xmax to x=0
						current = Coord(x_size - 1, center_y, (z_coeff * (z_size - 1)));
						target_x = 0;
						target_z = current.z;
						dx = -1;
						break;
						
					case 3: // z=zmax to z=0
						current = Coord((x_coeff * (x_size - 1)), center_y, z_size - 1);
						target_x = current.x;
						target_z = 0;
						dz = -1;
						break;
				}
				
				dy = 0;
				
				c_tumor.push_back(current);
				
				bool reached_end = false;
				int segment_length = 0;
				int max_length = 2 * std::max(x_size, z_size);  // Maximum path length
				
				while (!reached_end && segment_length < max_length) {
					// Gradually adjust direction with persistence
					double persistence = 0.2;
					if (rng.get_unif_01() > persistence) {
						if (seg % 4 == 0) {  // x=0 to x=max
							dx = 1;  // Keep moving right
							dz = (rng.get_unif_01() > 0.5) ? 1 : -1;  // Small z variation
						} 
						else if (seg % 4 == 1) {  // z=0 to z=max
							dx = (rng.get_unif_01() > 0.5) ? 1 : -1;  // Small x variation
							dz = 1;  // Keep moving forward
						}
						else if (seg % 4 == 2) {  // x=max to x=0
							dx = -1;  // Keep moving left
							dz = (rng.get_unif_01() > 0.5) ? 1 : -1;  // Small z variation
						}
						else {  // z=max to z=0
							dx = (rng.get_unif_01() > 0.5) ? 1 : -1;  // Small x variation
							dz = -1;  // Keep moving backward
						}
						
						// Y direction changes remain the same
						double r = rng.get_unif_01();
						if (r < 0.333) {
							dy = 1;
						} else if (r < 0.667) {
							dy = -1;
						} else {
							dy = 0;
						}
					}
					
					// Calculate next position
					int next_x = current.x + round(dx);
					int next_y = current.y + round(dy);
					int next_z = current.z + round(dz);
					
					// Bound checking and correction
					if (next_x < 0) {
						next_x = 0;
						dx = std::abs(dx);
					}
					if (next_x >= x_size) {
						next_x = x_size - 1;
						dx = -std::abs(dx);
					}
					
					if (next_y < 0) {
						next_y = 0;
						dy = 0;
					}
					if (next_y >= y_size) {
						next_y = y_size - 1;
						dy = 0;
					}
					
					// Special z-coordinate handling based on segment type
					if (seg % 4 == 0 || seg % 4 == 2) {  // Moving in x-direction
						// Keep z near its initial value (z_coeff * z_size)
						double target_z = z_coeff * z_size;
						if (std::abs(next_z - target_z) > 2) {  // Allow small variation
							next_z = current.z + (next_z > target_z ? -1 : 1);  // Move back toward target
						}
					} else {  // Moving in z-direction
						if (next_z < 0) {
							next_z = 0;
							dz = std::abs(dz);
						}
						if (next_z >= z_size) {
							next_z = z_size - 1;
							dz = -std::abs(dz);
						}
					}
					
					// Similarly for x-coordinate when moving in z-direction
					if (seg % 4 == 1 || seg % 4 == 3) {  // Moving in z-direction
						// Keep x near its initial value (x_coeff * x_size)
						double target_x = x_coeff * x_size;
						if (std::abs(next_x - target_x) > 2) {  // Allow small variation
							next_x = current.x + (next_x > target_x ? -1 : 1);  // Move back toward target
						}
					}
					
					// Update position
					current.x = next_x;
					current.y = next_y;
					current.z = next_z;
					
					// std::cout << "Current coord: " << current << ", dx: " << round(dx) << ", dy: " << round(dy) << ", dz: " << round(dz) << std::endl;
					// Check if we've reached target
					if ((seg % 4 == 0 && current.x >= target_x) ||
						(seg % 4 == 1 && current.z >= target_z) ||
						(seg % 4 == 2 && current.x <= target_x) ||
						(seg % 4 == 3 && current.z <= target_z)) {
						reached_end = true;
					}
					
					// Check distance from tumor center
					double d_x_center = current.x - center_x;
					double d_y_center = current.y - center_y;
					double d_z_center = current.z - center_z;
					double dist_sq = d_x_center * d_x_center + d_y_center * d_y_center + d_z_center * d_z_center;
					
					if (dist_sq > radius * radius) {
						c_tumor.push_back(current);
					}
					
					segment_length++;
				}
			}
		}
		for (size_t i = 0; i < c_tumor.size(); i++)
		{
			_tumor.add_init_vas(c_tumor[i]);
			//std::cout << "initial vas location: (" << c_tumor[i].x << ", " << c_tumor[i].y << ", " << c_tumor[i].z << ") " << std::endl;
			nr_source_tumor++;
		}
		
		std::cout << "Initial endothelial cells: " << c_tumor.size() << std::endl;
	}

	std::string s;
	_tumor.initCompartment(s);
	
}


/*! initialize compartments: create initial cells from file input specifications
	This function is called only when creating the model for the first time,
	and not reading from saved state. Objects created in this function should
	be already in the serialization and can be loaded directly.
*/
void HCC_Core::initializeSimulation(std::string core) {
	_tumor.set_allow_shift(false);
	_tumor.initCompartment(core);
}

void HCC_Core::preSimulation() {

	double t0 = 0;
	double tumor_volume = 0;
	const double sec_per_day = 86400;
	const double preSim_time = 5000;
	const double dt = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_SEC_PER_TIME_SLICE);
	std::cout << "Presimulation dt : " << dt << std::endl;
	const double tumor_volume_ref = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_INIT_TUM_VOL);
	double abm_min_cc = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_C1_MIN);

	const auto& stats = _tumor.get_stats();
	double init_cancer_count = double(stats.getCancerCell());
	std::cout << "populated cancer cells: " << init_cancer_count << std::endl;

	while (tumor_volume < tumor_volume_ref && t0 < preSim_time * sec_per_day) {
		std::cout << "Presimulation time: " << t0 / sec_per_day << std::endl;
		std::cout << "Presimulation tumor volume : " << tumor_volume << " , Reference Volume : " << tumor_volume_ref << std::endl;
		// std::cout << "RNG check (" << slice << ") START : " << rng.get_unif_01() << std::endl;

		/* update cancer number and blood concentration */
		auto& qsp_var = _lymph.get_var_exchange();
		double lymphCC = qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_C];

		// Length of the edge in the ABM module
		const double abm_length = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_VOXEL_SIZE_CM) * params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X) / 100;
		const double abm_volume = std::pow(abm_length, 3);

		/* if QSP halted, skip*/
		std::cout << "Presimulation lymph CC: " << lymphCC << std::endl;


		_tumor.update_abm_with_qsp(qsp_var);
		// Should equal to 0 at this point
		std::cout << "Pre treatment nivo: " << qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_NIVO] << std::endl;
		std::cout << "Pre treatment cabo: " << qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_CABO] << std::endl;
		std::cout << "Pre treatment ipi: " << qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_IPI] << std::endl;
		

		/* ABM time step */
		_tumor.timeSlice(t0);
		//std::cout << "RNG check (" << slice << ") MARGI : " << rng.get_unif_01() << std::endl;

		double qsp_tumor_volume = _tumor.get_Tum_Vol();

		/* update QSP variables */
		auto& abm_var_0 = _tumor.get_var_exchange();

		size_t abm_var_len = abm_var_0.size();
		auto abm_var = std::vector<double>(abm_var_len, 0);

		double w = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_WEIGHT_QSP);

		double tumCC = abm_var_0[SP_QSP_IO::SP_QSP_HCC::Tumor::TUMEX_CC];

		double abm_scaler = (1 - w) / w * lymphCC / (tumCC+ abm_min_cc );
		//double abm_scaler = (1 - w) / w * qsp_tumor_volume / abm_volume;

		std::cout << "Preliminary simulation scalor:\n" << abm_scaler << std::endl;
		for (size_t i = 0; i < abm_var_len; i++)
		{
			abm_var[i] = abm_var_0[i] * abm_scaler;
		}

		_lymph.update_qsp_var(abm_var);

		tumor_volume = _tumor.get_Tum_Vol();
		/* QSP time step */
		_lymph.time_step_preSimulation(t0, dt);

		t0 += dt;
		briefStats(t0);

		const auto& stats = _tumor.get_stats();
		int cancer_count_pretreatment = stats.getCancerCell();

		if (cancer_count_pretreatment < init_cancer_count / 2) {
			std::cout << "Cancer cell decreased abnormally without treatment." << std::endl;
			exit(0);
		}
	}
	if (tumor_volume < tumor_volume_ref)
	{
		std::cout << "tumor volume condition is not met" << std::endl;
		exit(0);
	}

	return;
}


void HCC_Core::timeSlice(const long slice){
	
	const double dt = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_SEC_PER_TIME_SLICE);
	const double t0 = slice * dt;
	std::cout << slice << std::endl;
	// std::cout << "RNG check (" << slice << ") START : " << rng.get_unif_01() << std::endl;

	/* update cancer number and blood concentration */
	auto& qsp_var = _lymph.get_var_exchange();
	double lymphCC = qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_C];

	// Length of the edge in the ABM module
	const double abm_length = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_VOXEL_SIZE_CM) * params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_TUMOR_X) / 100;
	const double abm_volume = std::pow(abm_length, 3);

	/* if QSP halted, skip*/
	std::cout << "lymph CC: " << lymphCC << std::endl;
	double abm_min_cc = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_C1_MIN);

	if (lymphCC > abm_min_cc)
	{
		_tumor.update_abm_with_qsp(qsp_var);
		std::cout << "nivo: " << qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_NIVO] << std::endl;
		std::cout << "cabo: " << qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_CABO] << std::endl;
		std::cout << "ipi: " << qsp_var[SP_QSP_IO::SP_QSP_HCC::LymphCentral::QSPEX_TUM_IPI] << std::endl;

		/*
		for (auto& v : qsp_var)
		{
		std::cout << v << ", ";
		}
		std::cout << std::endl;
		*/

		/* ABM time step */
		_tumor.timeSlice(slice);
		//std::cout << "RNG check (" << slice << ") MARGI : " << rng.get_unif_01() << std::endl;
		double qsp_tumor_volume = _tumor.get_Tum_Vol();
		/* update QSP variables */
		auto& abm_var_0 = _tumor.get_var_exchange();
		size_t abm_var_len = abm_var_0.size();
		auto abm_var = std::vector<double>(abm_var_len, 0);

		double w = params.getVal(SP_QSP_IO::SP_QSP_HCC::PARAM_WEIGHT_QSP);

		double tumCC= abm_var_0[SP_QSP_IO::SP_QSP_HCC::Tumor::TUMEX_CC];

		double abm_scaler = (1 - w) / w * lymphCC / (tumCC+ abm_min_cc );
		//double abm_scaler = (1 - w) / w * qsp_tumor_volume / abm_volume;

		std::cout << "scalor:\n" <<  abm_scaler<< std::endl;

		for (size_t i = 0; i < abm_var_len; i++)
		{
			abm_var[i] = abm_var_0[i] * abm_scaler;
		}

		_lymph.update_qsp_var(abm_var);

	}

	/* QSP time step */
	_lymph.time_step(t0, dt);
	
	return;

}

void HCC_Core::write_stats_header(void) const {

	auto& statsStream= output_hub.getStatsFstream();
	statsStream<< _tumor.get_stats().writeHeader();
	return;
}

void HCC_Core::write_stats_slice(unsigned long slice)const{
	{
		auto& statsStream = output_hub.getStatsFstream();
		statsStream << _tumor.get_stats().writeSlice(slice);
		statsStream.flush();
	}
	return;
}


void HCC_Core::write_QSP(unsigned long slice, bool header)const{
	auto& stream = output_hub.get_lymph_blood_QSP_stream();
	if (header){
		stream << "time" << _lymph.getQSPHeaders() << std::endl;
	}
	else{
		stream << slice << _lymph << std::endl;
	}
	return;
}


void HCC_Core::writeOde(unsigned long slice){
	_tumor.printCellOdeToFile(slice);
}


void HCC_Core::briefStats(unsigned long slice){
	std::cout << "Time: " << slice << std::endl;
	{
		const auto& stats = _tumor.get_stats();
		std::cout << "Core: " << "nrCell: " << _tumor.getNrCell()
			<< ", CD8: " << stats.getTCell()
			<< ", Th: " << stats.getTh()
			<< ", Treg: " << stats.getTreg()
			<< ", MDSC: " << stats.getMDSC()
			<< ", Cancer cell:" << stats.getCancerCell() 
			<< ", M1 Macrophage cell:" << stats.getMacM1()
			<< ", M2 Macrophage cell:" << stats.getMacM2()
			<< ", Vas cell:" << stats.getVas()
			<< ", APC: " << stats.getAPC()
			<< ", normal fibroblast: " << stats.getFibNormal()
			<< ", CAFs: " << stats.getFibCAFs() << std::endl;
	}
}


/*! Print grid info to file.
    \param [in] slice
	\param [in] option: 1. only cellular and vasculature scale; 2. only molecular scale; 3. both scales
*/
void HCC_Core::writeGrids(unsigned long slice, unsigned int option){
	if (option == 1)
	{
		{
			std::ofstream& cell = output_hub.getNewGridToSnapshotStream(slice, "cell_");
			cell << _tumor.compartment_cells_to_string();
			cell.close();

			std::ofstream& stroma = output_hub.getNewGridToSnapshotStream(slice, "stroma_");
			stroma << _tumor.printStromaToFile();
			stroma.close();
		}
		
	}
	if (option == 2)
	{
		{
			std::ofstream&  snap = output_hub.getNewGridToSnapshotStream(slice, "grid_core_");
			snap << _tumor.printGridToFile();
			snap.close();

			std::ofstream& grad_x = output_hub.getNewGridToSnapshotStream(slice, "gradient_x_");
			grad_x << _tumor.printGradientToFile(0);
			grad_x.close();

			std::ofstream& grad_y = output_hub.getNewGridToSnapshotStream(slice, "gradient_y_");
			grad_y << _tumor.printGradientToFile(1);
			grad_y.close();

			std::ofstream& grad_z = output_hub.getNewGridToSnapshotStream(slice, "gradient_z_");
			grad_z << _tumor.printGradientToFile(2);
			grad_z.close();

		}
	}


	if (option == 3)
	{
		{
			std::ofstream& cell = output_hub.getNewGridToSnapshotStream(slice, "cell_");
			cell << _tumor.compartment_cells_to_string();
			cell.close();

			std::ofstream& stroma = output_hub.getNewGridToSnapshotStream(slice, "stroma_");
			stroma << _tumor.printStromaToFile();
			stroma.close();

			std::ofstream& snap = output_hub.getNewGridToSnapshotStream(slice, "grid_core_");
			snap << _tumor.printGridToFile();
			snap.close();

		}
	}
}

