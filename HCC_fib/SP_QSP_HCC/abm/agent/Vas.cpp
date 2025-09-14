//#include <boost/serialization/export.hpp>
#include "Vas.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(Vas)
#include <iostream>
#include <sstream>

#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "../compartment/LymphCentral.h"

#include <iostream>
#include <sstream>

namespace SP_QSP_IO {
namespace SP_QSP_HCC {

	Vas::Vas(SpatialCompartment* c)
		:Cell_Tumor(c)
		, _source_O2(NULL)
		, _sink_VEGFA(NULL)
		, _divide_cd(0)
		, _tumble(true)
		, _moveDirection(0, 0, 0)
		, _upstream_neighbor(nullptr)
		, _tip_id(0)
		, _branch(false)
	{
		_state = AgentStateEnum::VAS_PHALANX;

	}

	Vas::Vas(const Vas& c)
		: Cell_Tumor(c)
		, _source_O2(NULL)
		, _sink_VEGFA(NULL)
		, _divide_cd(c._divide_cd)
	    , _divide_limit(c._divide_limit)
		, _tumble(true)
		, _moveDirection(c._moveDirection)
		, _upstream_neighbor(c._upstream_neighbor)
		, _tip_id(c._tip_id)
		, _branch(false)
	{
		setup_chem_sink(_sink_VEGFA, CHEM_VEGFA, params.getVal(PARAM_VEGFA_UPTAKE));
	}

	Vas::~Vas()
	{
	}

	std::string Vas::toString()const {
		std::stringstream ss;
		ss << Cell_Tumor::toString();
		return ss.str();
	}


	bool Vas::agent_movement_step(double t, double dt, Coord& c) {
		bool move = false;
		
		double ECM_density = get_tumor().get_coord_ECM(c);
		double ECM_sat = ECM_density / (ECM_density + params.getVal(PARAM_FIB_ECM_MOT_EC50));
		if (rng.get_unif_01() < ECM_sat){
			return move;
		}

		if (_state == AgentStateEnum::VAS_TIP)
		{
			//if (_divide_cd == 0) {
			if (_tumble == false) {
				//calculate the probability of tumble according to tgfb

				double v_x = _moveDirection.x / params.getVal(PARAM_SEC_PER_TIME_SLICE);
				double v_y = _moveDirection.y / params.getVal(PARAM_SEC_PER_TIME_SLICE);
				double v_z = _moveDirection.z / params.getVal(PARAM_SEC_PER_TIME_SLICE);


				std::vector< std::vector<double> > gradient = get_tumor().get_chem_grid().get_gradient(getCoord());

				/*!
					get gradient of(x,y,z); with given dimension d; x, y, z  are integer voxel coordinates
					d=0: x axis
					d=1: y axis
					d=2: z axis
					substrate i
				*/

				double dot_product = v_x * gradient[CHEM_VEGFA][0] + v_y * gradient[CHEM_VEGFA][1] + v_z * gradient[CHEM_VEGFA][2];
				//double delta = 5e-6;
				double lambda = params.getVal(PARAM_VAS_TUMBLE);

				double norm_v = std::sqrt(v_x * v_x + v_y * v_y + v_z * v_z);
				double norm_gradient = std::sqrt(gradient[CHEM_VEGFA][0] * gradient[CHEM_VEGFA][0] + gradient[CHEM_VEGFA][1] * gradient[CHEM_VEGFA][1] + gradient[CHEM_VEGFA][2] * gradient[CHEM_VEGFA][2]);
				if (norm_gradient == 0) {
					_tumble = true;
					return move;
				}
				double cos_theta = dot_product / (norm_v * norm_gradient);
				//std::cout << "cos_theta: " << cos_theta << ", nrom_v: " << norm_v << ", norm_gradient: " << norm_gradient << std::endl;
				double EC50_grad = 1;
				//Hill function of cell sensing the chemokine gradient, higher gradient makes cell easier to move (less likely to tumble)
				double H_grad = norm_gradient / (norm_gradient + EC50_grad);
				// if movement direction is opposite with
				if (cos_theta < 0) {
					H_grad = -H_grad;
				}
				double Tumble_rate = (lambda / 2) * (1 - cos_theta) * (1 - H_grad) * params.getVal(PARAM_SEC_PER_TIME_SLICE) + params.getVal(PARAM_VAS_DELTA); //0.2 is the random tumble coefficient/
				//double Tumble_rate = 0.3;
				// poisson process (tumble times N >= 1)
				double pTumble = 1 - std::exp(-Tumble_rate);
				//std::cout << "VEGF concentration: " << get_tumor().get_chem(getCoord(), CHEM_VEGFA) << ", x grad: " << gradient[CHEM_VEGFA][0] << ", y grad: " 
				//	<< gradient[CHEM_VEGFA][1] << ", z grad: " << gradient[CHEM_VEGFA][2] << std::endl;
				//std::cout << "v_x: " << v_x << ", v_y: " << v_y << ", v_z: " << v_z << ", dot_product: " << dot_product
				//	<< ", norm_v: " << norm_v << ", norm_gradient: " << norm_gradient << ", cos_theta: " << cos_theta
				//	<< ", EC50_grad: " << EC50_grad << ", H_grad: " << H_grad << ", Tumble_rate: " << Tumble_rate << ", pTumble: " << pTumble << std::endl;

				if (rng.get_unif_01() < pTumble)
				{
					_tumble = true;
					//std::cout << "cell tumbled, no move" << std::endl;
					return move;
				}

				Coord3D target_voxel = _moveDirection + getCoord();
				if (_compartment->voxelIsOpen(target_voxel, getType())) {
					move = true;
					c = target_voxel;
					//std::cout << "Run: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << target_voxel << std::endl;
				}
				return move;
			}
			else
			{
				std::vector<Coord3D> direction_vec;
				std::vector<double> prob_vec;
				double cdf = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{
						for (int k = -1; k < 2; k++)
						{
							if (i || j || k) {
								double dot_product = i * _moveDirection.x + j * _moveDirection.y + k * _moveDirection.z;
								double norm_dir = std::sqrt(i*i + j*j + k*k);
								double norm_movedir = std::sqrt(_moveDirection.x* _moveDirection.x + _moveDirection.y * _moveDirection.y + _moveDirection.z * _moveDirection.z);

								double cos_theta = dot_product / (norm_dir * norm_movedir);
								if (cos_theta > 0) {
									
									double sigma = 0.524;
									double rho = std::exp(cos_theta / sigma * sigma) / std::exp(1 / sigma * sigma);
									cdf += rho;
									direction_vec.push_back(Coord3D(i, j, k));
									prob_vec.push_back(cdf);
									//std::cout << "MoveDirection : " << _moveDirection << ", direction_vec: " << Coord3D(i, j, k) << ", dot_product: " << dot_product
									//	<< ", _moveDirection.x: " << _moveDirection.x << ", _moveDirection.y: " << _moveDirection.y << ", _moveDirection.z: " << _moveDirection.z
									//	<< ", cos_theta: " << cos_theta << std::endl;
								}
							}
						}
					}
				}
				for (double& element : prob_vec) {
					element /= cdf; // Divide each element by the divisor
				}

				//std::cout << "Current position: " << getCoord() << ",Current direction: " << _moveDirection << ", Candidate Direction: " 
				//<< ",length of prob_vec: " << prob_vec.size() << std::endl;
				
				for (int i = 0; i < prob_vec.size(); i++) {
					//std::cout << direction_vec[i] << ", Prob: " << prob_vec[i];
				}
				// Safety check - ensure we have valid directions
				if (direction_vec.empty() || prob_vec.empty()) {
					_tumble = false;  // Reset tumble state
					move = false;
					return move;     // No movement this step
				}
				size_t idx = rng.sample_cdf(prob_vec);
				Coord3D target_voxel = direction_vec[idx] + getCoord();
				
				if (_compartment->voxelIsOpen(target_voxel, getType())) {
					move = true;
					c = target_voxel;
					_moveDirection = direction_vec[idx];
					_tumble = false;
					//std::cout << "Run: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << target_voxel << std::endl;
				}

				// move
				/*
				int idx;
				const auto shape = getCellShape();
				if (_compartment->getOneOpenVoxel(shape->getVasDestinationVoxels(),
					shape->getVasDirectionAnchor(), getCoord(), getType(), idx, rng))
				{
					move = true;
					c = getCellShape()->getVasDirectionAnchor()[idx] + getCoord();
					_moveDirection = getCellShape()->getVasDirectionAnchor()[idx];
					//std::cout << "Tumble: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << c << std::endl;
					_tumble = false;
				}
				*/
			}
			//}
		}
		return move;
	}

	void Vas::make_mature() {
		Vas* v = this;
		while (v) {
			Vas* prev = v;
			v->set_phalanx();
			v = v->get_upstream_neighbor();
			prev->set_upstream_neighbor(nullptr);
		}
	}

	bool Vas::agent_state_step(double t, double dt, Coord& c) {
		const auto shape = getCellShape();
		bool divide = false;
		
		if (_state == AgentStateEnum::VAS_TIP)
		{
			//std::cout << "Tip move direction: " << _moveDirection << std::endl;
			/*
			_compartment->for_each_neighbor_ag(shape->getEnvironmentLocations(),
				getCoord(), [&](BaseAgent* ag) {
					auto pCell = dynamic_cast<Cell_Tumor*>(ag);
					if (ag->getType() == AgentTypeEnum::CELL_TYPE_VAS) {
						Vas* pVas_neighbor = dynamic_cast<Vas*>(ag);
						if ((pVas_neighbor->getState() == AgentStateEnum::VAS_STALK || pVas_neighbor->getState() == AgentStateEnum::VAS_TIP) && get_tip_id() != pVas_neighbor->get_tip_id()) {
							// we met another stalk, time to make them mature
							pVas_neighbor->make_mature();
							make_mature();
						}
					}
					return false;
				});
		    
			std::cout << "Tip cell at: " << getCoord() << " divide" << std::endl;
			*/
			//Coord3D target_voxel = _moveDirection + getCoord();
			if (_compartment->voxelIsOpen(getCoord(), getType())) {
				//std::cout << "check vas availablility :  success: " << getCoord() << std::endl;
				c = getCoord();
				divide = true;
				//std::cout << "Run: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << target_voxel << std::endl;
			}
			else {
				//std::cout << "check vas availablility :  fail: " << getCoord() << std::endl;
			}
			return divide;
		}


		if (_state == AgentStateEnum::VAS_STALK)
		{
			
		}

		if (_state == AgentStateEnum::VAS_PHALANX)
		{
			
			double vegf = get_tumor().get_chem(getCoord(), CHEM_VEGFA);
			double p_tip_cell = vegf / (params.getVal(PARAM_VAS_50) + vegf);
			
			
			//std::cout << "VEGF: " << vegf << ", p tip cell: " << p_tip_cell << std::endl;
			
			if (rng.get_unif_01() < p_tip_cell) {
				CoordVec tip_neighbor_CoordVec;
				int num_neighbors = params.getVal(PARAM_VAS_MIN_NEIGHBOR);
				//std::cout << "generating coords for " <<getCoord() << " : ";
				for (int i = -num_neighbors; i < num_neighbors; i++) {
					for (int j = -num_neighbors; j < num_neighbors; j++) {
						for (int k = -num_neighbors; k < num_neighbors; k++) {
							if (i || j || k) {
								Coord c_neighbor(i, j, k);
								tip_neighbor_CoordVec.push_back(c_neighbor);
								//std::cout << "c_neighbor: "<< c_neighbor << ", neighbor coords: " << neighbor_coord << ". ";
							}
						}
					}
				}
				//std::cout << std::endl;
				double nearby_vas_cell_exist = false;
				//std::cout << "check nearby Vas cell for differntiation " << getCoord() <<". ";
				_compartment->for_each_neighbor_ag(tip_neighbor_CoordVec,
					getCoord(), [&](BaseAgent* ag) { 
						auto pCell = dynamic_cast<Cell_Tumor*>(ag);
						//std::cout << " nearby vas cell: " << pCell->getCoord();
						if (pCell->getType() == AgentTypeEnum::CELL_TYPE_VAS) {
							auto pVas = dynamic_cast<Vas*>(pCell);
							if (_tip_id != pVas->get_tip_id()) {
								nearby_vas_cell_exist = true;
								//std::cout << ", with different tip ID. " ;
							}
							else {
								//std::cout << ", with same tip ID. ";
							}
						}
						return true;
					});
				//std::cout << std::endl;
				//std::cout << "Vas cell at: " << getCoord() << " cleared" << std::endl;
				if (!nearby_vas_cell_exist) {
					c = getCoord();
					divide = true;
					_branch = false;
				}
				
			}
			
		}

		return divide;

	}
	
	void Vas::set_tip() {
		_state = AgentStateEnum::VAS_TIP;
	}

	void Vas::set_stalk() {
		_state = AgentStateEnum::VAS_STALK;
	}

	void Vas::set_phalanx() {
		_state = AgentStateEnum::VAS_PHALANX;

		double vegf = get_tumor().get_chem(getCoord(), CHEM_VEGFA);
		double p_branch = params.getVal(PARAM_VAS_BRANCH_PROB) * vegf / (params.getVal(PARAM_VAS_50) + vegf);
		if (rng.get_unif_01() < p_branch) {
			//std::cout << "cell: " << getCoord() << " trying to branch" << std::endl;
			set_branch();
		}

	}

	void Vas::set_upstream_neighbor(Vas* upstream_neighbor) {
		_upstream_neighbor = upstream_neighbor;
	}

	void Vas::set_tip_id(unsigned int tip_id) {
		_tip_id = tip_id;
	}

	void Vas::set_branch() {
		_branch = true;
	}
	void Vas::set_moveDirection() {
		int idx;
		const auto shape = getCellShape();
		if (_compartment->getOneOpenVoxel(shape->getVasDestinationVoxels(),
			shape->getVasDirectionAnchor(), getCoord(), getType(), idx, rng))
		{
			_moveDirection = getCellShape()->getVasDirectionAnchor()[idx];
		}
		Coord3D current_position = getCoord();
		if (_moveDirection.x == 0 && _moveDirection.y == 0 && _moveDirection.z == 0) {
            _moveDirection.x = 1;  // Default movement if all zeros
        }
		if (current_position.x < 0.1 * params.getVal(PARAM_TUMOR_X)) {
			_moveDirection.x = 1;
		}
		if (current_position.x > 0.9 * params.getVal(PARAM_TUMOR_X)) {
			_moveDirection.x = -1;
		}
		if (current_position.y < 0.1 * params.getVal(PARAM_TUMOR_Y)) {
			_moveDirection.y = 1;
		}
		if (current_position.y > 0.9 * params.getVal(PARAM_TUMOR_Y)) {
			_moveDirection.y = -1;
		}
		if (current_position.z < 0.1 * params.getVal(PARAM_TUMOR_Z)) {
			_moveDirection.z = 1;
		}
		if (current_position.z > 0.9 * params.getVal(PARAM_TUMOR_Z)) {
			_moveDirection.z = -1;
		}
		
		
		//std::cout << "set move direction: " << _moveDirection << std::endl;
	}
	void Vas::remove_all_source_sink(void) {

		remove_source_sink(_source_O2);
		remove_source_sink(_sink_VEGFA);

		return;
	}
};
};