//#include <boost/serialization/export.hpp>
#include "Mac.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(Mac)
#include <iostream>
#include <sstream>

#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "../compartment/LymphCentral.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

Mac::Mac(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _source_IFNg(NULL)
	, _source_IL_12(NULL)
	, _source_TGFB(NULL)
	, _source_IL_10(NULL)
	, _source_VEGFA(NULL)
	, _sink_CCL2(NULL)
{
	_state = AgentStateEnum::MAC_M1;
	_life = getMacLife();
	
}

Mac::Mac(const Mac& c)
	:Cell_Tumor(c)
	, _source_IFNg(NULL)
	, _source_IL_12(NULL)
	, _source_TGFB(NULL)
	, _source_IL_10(NULL)
	, _source_VEGFA(NULL)
	, _sink_CCL2(NULL)
{
	_life = getMacLife();
	if (_state == AgentStateEnum::MAC_M1)
	{
		//once M1 macrophage in contact with cancer cell, they release proinflammatory cytokine
		setup_chem_source(_source_IFNg, CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
		setup_chem_source(_source_IL_12, CHEM_IL_12, params.getVal(PARAM_IL_12_RELEASE));
		setup_chem_source(_source_TGFB, CHEM_TGFB, 0);
		setup_chem_source(_source_IL_10, CHEM_IL_10, 0);
		setup_chem_source(_source_VEGFA, CHEM_VEGFA, 0);
		//sink
		setup_chem_sink(_sink_CCL2, CHEM_CCL2, params.getVal(PARAM_CCL2_UPTAKE));

	}
	if (_state == AgentStateEnum::MAC_M2)
	{
		setup_chem_source(_source_IFNg, CHEM_IFN, 0);
		setup_chem_source(_source_IL_12, CHEM_IL_12, 0);
		setup_chem_source(_source_TGFB, CHEM_TGFB, params.getVal(PARAM_MAC_TGFB_RELEASE));
		setup_chem_source(_source_IL_10, CHEM_IL_10, params.getVal(PARAM_MAC_IL_10_RELEASE));
		setup_chem_source(_source_VEGFA, CHEM_VEGFA, params.getVal(PARAM_MAC_VEGFA_RELEASE));
		//sink
		setup_chem_sink(_sink_CCL2, CHEM_CCL2, params.getVal(PARAM_CCL2_UPTAKE));
	}

}

Mac::~Mac()
{
}

std::string Mac::toString()const {
	std::stringstream ss;
	ss << Cell_Tumor::toString();
	return ss.str();
}

bool Mac::agent_movement_step(double t, double dt, Coord& c) {

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

bool Mac::agent_state_step(double t, double dt, Coord& c) {
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

	if (!isDead() && _state == AgentStateEnum::MAC_M1)
	{

		double TGFB = get_tumor().get_chem(getCoord(), CHEM_TGFB);
		double IL10 = get_tumor().get_chem(getCoord(), CHEM_IL_10);

		//! calculate M2 polarization constant alpha based on TGFB and IL10 secretion
		double alpha = params.getVal(PARAM_MAC_M2_POL) * ((TGFB / (TGFB + params.getVal(PARAM_MAC_TGFB_EC50))) + (IL10 / (IL10 + params.getVal(PARAM_MAC_IL_10_EC50))));
		//double alpha = params.getVal(PARAM_MAC_M2_POL);
		//! calculate probability of M2 polarization based on alpha
		double p_M2_polar = 1 - std::exp(-alpha);

		/*
		std::cout 
			<< "Coord: " << getCoord() << ", M1 to M2 Polarization constant: " << params.getVal(PARAM_MAC_M2_POL)
			<< ", TGFB : " << TGFB << ", PARAM_MAC_TGFB_EC50: " << params.getVal(PARAM_MAC_TGFB_EC50) << ", H_TGFB: " << (TGFB / (TGFB + params.getVal(PARAM_MAC_TGFB_EC50)))
			<< ", IL10: " << IL10 << ", PARAM_MAC_IL_10_EC50: " << params.getVal(PARAM_MAC_IL_10_EC50) << ", H_IL10: " << (IL10 / (IL10 + params.getVal(PARAM_MAC_IL_10_EC50)))
			<< ", M1 to M2 Polarization constant alpha: " << alpha << ", probability p_M2_polar: " << p_M2_polar << std::endl;
		*/

		if (rng.get_unif_01() < p_M2_polar)
		{
			setM2();
			//std::cout << "M1 to M2 macrophage polarization" << std::endl;
			return divide;
		}

	}
	if (!isDead() && _state == AgentStateEnum::MAC_M2)
	{
		double IL12 = get_tumor().get_chem(getCoord(), CHEM_IL_12);
		double IFNg = get_tumor().get_chem(getCoord(), CHEM_IFN);

		//! calculate M1 polarization constant alpha based on IL12 and IFNg secretion
		double alpha = params.getVal(PARAM_MAC_M1_POL) * ((IFNg / (IFNg + params.getVal(PARAM_MAC_IFN_G_EC50))) + (IL12 / (IL12 + params.getVal(PARAM_MAC_IL_12_EC50))));
		//! calculate probability of M1 polarization based on alpha
		double p_M1_polar = 1 - std::exp(-alpha);

		/*
		std::cout 
			<< "M2 to M1 Polarization constant: " << params.getVal(PARAM_MAC_M1_POL)
			<< ", IFNg : " << IFNg << ", PARAM_MAC_IFN_G_EC50: " << params.getVal(PARAM_MAC_IFN_G_EC50) << ", H_IFNg: " << (IFNg / (IFNg + params.getVal(PARAM_MAC_IFN_G_EC50)))
			<< ", IL12: " << IL12 << ", PARAM_MAC_IL_12_EC50: " << params.getVal(PARAM_MAC_IL_12_EC50) << ", H_IL12: " << (IL12 / (IL12 + params.getVal(PARAM_MAC_IL_12_EC50)))
			<< ", M2 to M1 Polarization constant alpha: " << alpha << ", probability p_M1_polar: " << p_M1_polar << std::endl;
		*/

		if (rng.get_unif_01() < p_M1_polar)
		{
			setM1();
			//std::cout << "M2 to M1 macrophage polarization" << std::endl;
			return divide;
		}
	}
	return divide;
}

void Mac::agent_state_scan(void) {

	const auto shape = getCellShape();
	if (_state == AgentStateEnum::MAC_M1) {
		/**/
		// scan cancer cells
		//std::cout << "T cell at: "<< getCoord() << std::endl;
		//int nr_cancer_neighbor = 0;
		_compartment->for_each_neighbor_ag(shape->getEnvironmentLocations(),
			getCoord(), [&](BaseAgent* ag) {
				//_count_neighbor_all++;
				auto pCell = dynamic_cast<Cell_Tumor*>(ag);
				//update_neighbor_PDL1(pCell->get_PDL1());
				if (ag->getType() == AgentTypeEnum::CELL_TYPE_CANCER) {
					auto pCancer = dynamic_cast<CancerCell*>(ag);
					pCancer->inc_neighbor_MacM1();
					//once M1 macrophage in contact with cancer cell, they release proinflammatory cytokine
					update_chem_source(_source_IFNg, params.getVal(PARAM_IFN_G_RELEASE));
					update_chem_source(_source_IL_12, params.getVal(PARAM_IL_12_RELEASE));
					//nr_cancer_neighbor += 1;
				}
				return false;
			});
		//std::cout << "T cell neighbor PDL1: "<< _max_neighbor_PDL1 << std::endl;
		//std::cout << "T cell neighbor Cancer: "<< nr_cancer_neighbor << std::endl;
	}
	return;
}

double Mac::get_PD1_PDL1(double PDL1, double Nivo) {

	static double T1 = params.getVal(PARAM_MAC_PD1_SYN);
	double T2 = PDL1;
	static double k1 = params.getVal(PARAM_PDL1_K1);
	static double k2 = params.getVal(PARAM_PDL1_K2);
	static double k3 = params.getVal(PARAM_PDL1_K3);
	/*
	double c = T1 + PDL1 + (1 + k2 * Nivo) / k1;
	double discriminant = c * c - 4 * T1 * PDL1;
	double bond = .5*(c - std::sqrt(discriminant));
	*/

	double a = 1;
	double b = (Nivo * k2 / k1 * (2 * k3 / k1 - 1) - 2 * T2 - T1 - 1 / k1) / T2;
	double c = (Nivo * k2 / k1 + 1 / k1 + T2 + 2 * T1) / T2;
	double d = -T1 / T2;

	/*
	std::cout << k1 << "," << k2 << "," << k3 << ","
		<< T1 << "," << T2 << std::endl;
	std::cout << a << "," << b << "," << c << "," << d << std::endl;
	*/

	//Newton_Raphson_root
	int max_iter = 20;
	double tol_rel = 1E-5;
	double root = 0;
	double res, root_new, f, f1;
	int i = 0;
	while (i < max_iter) {
		f = a * std::pow(root, 3) + b * std::pow(root, 2) + c * root + d;
		f1 = 3.0 * a * std::pow(root, 2) + 2.0 * b * root + c;
		root_new = root - f / f1;
		res = std::abs(root_new - root) / root_new;
		if (res > tol_rel) {
			i++;
			root = root_new;
		}
		else {
			break;
		}
	}

	
	//std::cout << "PDL1: " << PDL1 << ", nivo: " << Nivo << ", root: " << root << std::endl;
	

	double bond = T2 * root;
	//std::cout << "bond : " << bond << std::endl;

	return bond;
}
//get Hill function result for binding molecules
double Mac::get_H(double bond, double n, double k50) {
	return get_Hill_equation(bond, k50, n);
}


double Mac::get_kill_prob(double p, double q) {
	return 1 - std::pow(p, q);
}



void Mac::move_all_source_sink(void) const {
	move_source_sink(_source_IFNg);
	move_source_sink(_source_IL_12);
	move_source_sink(_source_TGFB);
	move_source_sink(_source_IL_10);
	move_source_sink(_source_VEGFA);
	move_source_sink(_sink_CCL2);
}

void Mac::remove_all_source_sink(void) {
	remove_source_sink(_source_IFNg);
	remove_source_sink(_source_IL_12);
	remove_source_sink(_source_TGFB);
	remove_source_sink(_source_IL_10);
	remove_source_sink(_source_VEGFA);
	move_source_sink(_sink_CCL2);
}

void Mac::setM1() {
	_state = AgentStateEnum::MAC_M1;

	//! Change cytokine secretion
	update_chem_source(_source_IFNg, params.getVal(PARAM_IFN_G_RELEASE));
	update_chem_source(_source_IL_12, params.getVal(PARAM_IL_12_RELEASE));
	update_chem_source(_source_TGFB, 0.0);
	update_chem_source(_source_IL_10, 0.0);
	update_chem_source(_source_VEGFA, 0.0);

}

void Mac::setM2() {
	_state = AgentStateEnum::MAC_M2;
	
	//! Change cytokine secretion
	update_chem_source(_source_TGFB, params.getVal(PARAM_MAC_TGFB_RELEASE));
	update_chem_source(_source_IL_10, params.getVal(PARAM_MAC_IL_10_RELEASE));
	update_chem_source(_source_VEGFA, params.getVal(PARAM_MAC_VEGFA_RELEASE));
	update_chem_source(_source_IFNg, 0);
	update_chem_source(_source_IL_12, 0);
}
int Mac::getMacLife() {
	double lifeMean = params.getVal(PARAM_MAC_LIFE_MEAN);

	double tLifeD = rng.get_exponential(lifeMean);

	int tLife = int(tLifeD + 0.5);
	tLife = tLife > 0 ? tLife : 0;
	//std::cout << "random TCD4 life: " << tLife << std::endl;
	return tLife;
}

};
};