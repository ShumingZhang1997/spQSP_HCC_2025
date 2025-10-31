//#include <boost/serialization/export.hpp>
#include "MDSC.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(MDSC)

#include <iostream>
#include <sstream>

#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "../compartment/LymphCentral.h"
//#include "TCell.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

MDSC::MDSC(SpatialCompartment* c)
	:Cell_Tumor(c)
	//, _source_IL_10(NULL)
	, _source_ArgI(NULL)
	, _source_NO(NULL)    
{
	_life = getMDSCLife();
}

MDSC::MDSC(const MDSC& c)
	:Cell_Tumor(c)
	//, _source_IL_10(NULL)
	, _source_ArgI(NULL)
	, _source_NO(NULL)     
	, _sink_CCL2(NULL)
{
	_life = getMDSCLife();
	//source
	setup_chem_source(_source_ArgI, CHEM_ARGI, params.getVal(PARAM_ARGI_RELEASE));
	setup_chem_source(_source_NO, CHEM_NO, params.getVal(PARAM_NO_RELEASE));
	//sink
	setup_chem_sink(_sink_CCL2, CHEM_CCL2, params.getVal(PARAM_CCL2_UPTAKE));
}

MDSC::~MDSC()
{
}

std::string MDSC::toString()const{
	std::stringstream ss;
	ss << Cell_Tumor::toString();
	return ss.str();
}

//void MDSC::setDead(void)
//{   
//	//remove_source_sink(_source_IL_10);
//}

bool MDSC::agent_movement_step(double t, double dt, Coord& c){
	
	bool move = false;
	double ECM_density = get_tumor().get_coord_ECM(c);
	double ECM_sat = ECM_density / (ECM_density + params.getVal(PARAM_FIB_ECM_MOT_EC50));
	if (rng.get_unif_01() < ECM_sat){
		return move;
	}
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

		double dot_product = v_x * gradient[CHEM_CCL2][0] + v_y * gradient[CHEM_CCL2][1] + v_z * gradient[CHEM_CCL2][2];
		
		//double delta = 5e-6;
		double lambda = 0.0000168;
	
		double norm_v = std::sqrt(v_x * v_x + v_y * v_y + v_z * v_z);
		double norm_gradient = std::sqrt(gradient[CHEM_CCL2][0] * gradient[CHEM_CCL2][0] + gradient[CHEM_CCL2][1] * gradient[CHEM_CCL2][1] + gradient[CHEM_CCL2][2] * gradient[CHEM_CCL2][2]);
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
	
	return move;
}

bool MDSC::agent_state_step(double t, double dt, Coord& c){
	bool divide = false;
	if (!isDead())
	{
		_life--;
		if (_life <= 0)
		{
			setDead();
			// remove source when cell die
			return divide;
		}
	}

	const auto shape = getCellShape();
	Cell_Tumor::agent_state_step(t, dt, c);

	return divide;
}


int MDSC::getMDSCLife(){

	double lifeMean = params.getVal(PARAM_MDSC_LIFE_MEAN);

	double tLifeD = rng.get_exponential(lifeMean);

	int tLife = int(tLifeD + 0.5);
	tLife = tLife > 0 ? tLife : 0;
	
	//std::cout << "PARAM_MDSC_LIFE_MEAN: " << lifeMean
		//<< ", random MDSC life: " << tLife << std::endl;
	
	return tLife;
}

//! move sources (NO and ArgI)
void MDSC::move_all_source_sink(void)const{
	//std::cout << "moving sources: " << getCoord() << std::endl;
	move_source_sink(_source_ArgI);
	move_source_sink(_source_NO);  
	move_source_sink(_sink_CCL2);
	return;
}

//! remove sources (NO and ArgI)
void MDSC::remove_all_source_sink(void){
	remove_source_sink(_source_ArgI);
	remove_source_sink(_source_NO); 
	move_source_sink(_sink_CCL2);
	return;
}

};
};