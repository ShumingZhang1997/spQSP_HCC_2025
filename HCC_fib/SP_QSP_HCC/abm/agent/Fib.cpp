//#include <boost/serialization/export.hpp>
#include "Fib.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(Fib)
#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "../compartment/Tumor.h"
#include "../../core/GlobalUtilities.h"

#include <iostream>
#include <sstream>

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

Fib::Fib(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _next(nullptr)
	, _previous(nullptr)
{
	_state = AgentStateEnum::FIB_NORMAL;
}

Fib::Fib(const Fib& c)
	:Cell_Tumor(c)
	, _next(nullptr)
	, _previous(nullptr)
{

}
Fib::~Fib()
{
}

bool Fib::agent_movement_step(double t, double dt, Coord& c) {

	bool move = false;
	double ECM_density = get_tumor().get_coord_ECM(c);
	double ECM_sat = ECM_density / (ECM_density + params.getVal(PARAM_FIB_ECM_MOT_EC50));
	if (rng.get_unif_01() < ECM_sat){
		return move;
	}
	// for fibroblasts, only the head should "move" here; the tail is dragged along with it by special-case code in Tumor.cpp
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

		double dot_product = v_x * gradient[CHEM_TGFB][0] + v_y * gradient[CHEM_TGFB][1] + v_z * gradient[CHEM_TGFB][2];
		
		//double delta = 5e-6;
		double lambda = 0.0000168;
		if (_state == AgentStateEnum::FIB_CAF) {
			double lambda = 0.000168;
		}
		double norm_v = std::sqrt(v_x * v_x + v_y * v_y + v_z * v_z);
		double norm_gradient = std::sqrt(gradient[CHEM_TGFB][0] * gradient[CHEM_TGFB][0] + gradient[CHEM_TGFB][1] * gradient[CHEM_TGFB][1] + gradient[CHEM_TGFB][2] * gradient[CHEM_TGFB][2]);
		double cos_theta = dot_product / (norm_v * norm_gradient);

		double EC50_grad = 1;
		//Hill function of cell sensing the chemokine gradient, higher gradient makes cell easier to move (less likely to tumble)
		double H_grad = norm_gradient / (norm_gradient + EC50_grad);
		// if movement direction is opposite with tgfb gradient
		if (cos_theta < 0) {
			H_grad = -H_grad;
		}
		double Tumble_rate = (lambda/2) * (1 - cos_theta) * (1 - H_grad) * params.getVal(PARAM_SEC_PER_TIME_SLICE);
		//double Tumble_rate = 0.3;
		// poisson process (tumble times N >= 1)
		double pTumble = 1 - std::exp(-Tumble_rate);
		//std::cout << "TGFB concentration: " << get_tumor().get_chem(getCoord(), CHEM_TGFB) << ", x grad: " << gradient[CHEM_TGFB][0] << ", y grad: " << gradient[CHEM_TGFB][1] << ", z grad: " << gradient[CHEM_TGFB][2] << std::endl;
		//std::cout << "v_x: " << v_x << ", v_y: " << v_y << ", v_z: " << v_z << ", dot_product: " << dot_product 
		//	<< ", norm_v: " << norm_v << ", norm_gradient: " << norm_gradient << ", cos_theta: " << cos_theta 
		//	<< ", EC50_grad: " << EC50_grad << ", H_grad: " << H_grad << ", Tumble_rate: " << Tumble_rate << ", pTumble: " << pTumble << std::endl;

		if (rng.get_unif_01() < pTumble)
		{
			_tumble = true;
			//std::cout << "cell tumbled, no move" << std::endl;
			return move;
		}
	}
	
	if (_previous == nullptr)
	{
		// when the cell decides if it moves it tracks the tumble status from last step
		if (_tumble == true)
		{
			// move
			int idx;
			const auto shape = getCellShape();
			if (_compartment->getOneOpenVoxel(shape->getMoveDestinationVoxels(),
				shape->getMoveDirectionAnchor(), getCoord(), getType(), idx, rng))
			{
				move = true;
				c = getCellShape()->getMoveDirectionAnchor()[idx] + getCoord();
				_moveDirection = getCellShape()->getMoveDirectionAnchor()[idx];
				//std::cout << "Tumble: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << c << std::endl;
				_tumble = false;
			}
		}
		else {
			Coord3D target_voxel = _moveDirection + getCoord();
			if (_compartment->voxelIsOpen(target_voxel, getType())) {
				move = true;
				c = target_voxel;
				//std::cout << "Run: Previous location: " << getCoord() << ", move direction: " << _moveDirection << ", new location: " << target_voxel << std::endl;
			}
		}
	}
	return move;
}

bool Fib::agent_state_step(double t, double dt, Coord& c) {
	bool expand = false;
	
	/*
	  fibroblast only expand upon activation when fibroblast turns to CAFs depends on TGFb concentation
	  it does not proliferate.
	  The fibroblast only expand at the last agent
	 */
	double TGFB = get_tumor().get_chem(getCoord(), CHEM_TGFB);
	double activation_coeff = params.getVal(PARAM_FIB_CAF_ACTIVATION) * 5 * (1 + (TGFB / (TGFB + params.getVal(PARAM_FIB_CAF_EC50))));
	//double activation_coeff = 0.05 * (1 + (TGFB / (TGFB + params.getVal(PARAM_FIB_CAF_EC50))));
	double p_activation = 1 - std::exp(-activation_coeff);
	//std::cout << "TGFB: " << TGFB << ", activation_coeff: " << activation_coeff << ", PARAM_FIB_CAF_ACTIVATION: " << params.getVal(PARAM_FIB_CAF_ACTIVATION) << ", p_activation: " << p_activation << std::endl;
	if (_next == nullptr && rng.get_unif_01() < p_activation && _state == AgentStateEnum::FIB_NORMAL)
	{
		expand = true;
	}
	return expand;
};

void Fib::setNext(Fib* next) {
	_next = next;
};

void Fib::setPrevious(Fib* previous) {
	_previous = previous;
};

void Fib::setCAF() {
	_state = AgentStateEnum::FIB_CAF;
};

//! move sources (CCL2)
void Fib::move_all_source_sink(void)const {
	//std::cout << "moving sources: " << getCoord() << std::endl;
	//move_source_sink(_source_FIB_FIELD);
	return;
}

//! remove sources (CCL2)
void Fib::remove_all_source_sink(void) {
	return;
}

std::string Fib::toString()const{
	std::stringstream ss;
	ss << Cell_Tumor::toString();
	return ss.str();
}
};
};

/*
// assume `this` is the tail
if (grow && _next == nullptr && _state == ) {
	if (Coord3D c1 = get_tumor().getOpenNeighborForFib(c)) {
		if (Coord3D c2 = get_tumor().getOpenNeighborForFib(c1)) {
			CellAgent* f1 = createOneInitCell(AgentTypeEnum::CELL_TYPE_FIB, state, c1);
			CellAgent* f2 = createOneInitCell(AgentTypeEnum::CELL_TYPE_FIB, state, c2);

			this->setNext(f1);
			f1->setPrevious(this);
			f1->setNext(f2);
			f2->setPrevious(f1);
		}
	}
}
*/