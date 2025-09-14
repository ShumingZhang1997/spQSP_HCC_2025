#include "VoxelContentGen.h"
#include "../../core/GlobalUtilities.h"

namespace SP_QSP_IO{
namespace SP_QSP_HCC{

VoxelContentGen::VoxelContentGen()
	: _stationary(true)
	, _x_min(0)
	, _y_min(0)
	, _z_min(0)
	, _x_max(0)
	, _y_max(0)
	, _z_max(0)
	, _celltype_cdf()
	, _radius(0)
	, _center_x(0)
	, _center_y(0)
	, _center_z(0)
	, _x_fib(-1)
	, _y_fib(-1)
	, _z_fib(-1)
	, _include_fib(false)
{
	int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
	// 0: stem; 1-dmax: progenitor; dmax+1: senescent; dmax+2: empty
	_celltype_cdf = std::vector<double>(dmax + 3, 0.0);
	return;
}


VoxelContentGen::~VoxelContentGen()
{
}

void VoxelContentGen::set_cell_cdf(bool stationary) {
	unsigned int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
	double k, r, rs, rp, mu, l0, l1, l2;
	k = params.getVal(PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB);
	rs = params.getVal(PARAM_CSC_GROWTH_RATE);
	rp = params.getVal(PARAM_CANCER_PROG_GROWTH_RATE);
	mu = params.getVal(PARAM_CANCER_SENESCENT_DEATH_RATE);
	r = rs * (1 - k);
	l0 = k*rs / (r + rp);
	l1 = 2 * rp / (r + rp);
	l2 = 2 * rp / (r + mu);
	double C;
	if (stationary) {
		if (l1 == 1) {
			C = 1 + l0 + l0 * l2 * std::pow(l1, (dmax - 1));
		}
		else {
			C = 1 + l0 * (std::pow(l1, dmax) - 1) / (l1 - 1) + l0 * l2 * std::pow(l1, (dmax - 1));
		}
		double p;
		_celltype_cdf[0] = p = 1 / C; // joint P
		p *= l0;
		_celltype_cdf[1] = _celltype_cdf[0] + p;

		for (size_t i = 2; i <= dmax; i++)
		{
			p *= l1;
			_celltype_cdf[i] = _celltype_cdf[i - 1] + p;
		}
		p *= l2;
		_celltype_cdf[dmax + 1] = _celltype_cdf[dmax] + p;
		_celltype_cdf[dmax + 2] = 1.0;
	}
	else {
		double pstem = params.getVal(PARAM_DENSITY_CSC);
		double p = (1 - pstem) / double(dmax + 1);
		std::cout << "p step: " << p << std::endl;
		_celltype_cdf[0] = pstem;

		for (size_t i = 1; i < dmax + 2; i++)
		{
			_celltype_cdf[i] = _celltype_cdf[i - 1] + p;
		}

		_celltype_cdf[dmax + 2] = 1.0;
	}
	/*
	std::cout << "cancer cell cdf:" << std::endl;
	for (size_t i = 0; i < _celltype_cdf.size(); i++)
	{
		std::cout << i << ", " << _celltype_cdf[i] << std::endl;
	}
	*/
	return;
}
/*! return true if any cell is to be created
*/
bool VoxelContentGen::get_type_state(const Coord3D& c, RNG& rng,
	AgentType& type, AgentState& state, int&div)const{

	bool create_cell = false;
	bool cell_type_sample = false;

	if (!_circular)
	{
		if (_stationary || ((c.x >= _x_min && c.x < _x_max)
			&& (c.y >= _y_min && c.y < _y_max)
			&& (c.z >= _z_min && c.z < _z_max)))
		{
			cell_type_sample = true;
		}
	}

	if (_circular) {
		double distance_square = (c.x - _center_x) * (c.x - _center_x) + (c.y - _center_y) * (c.y - _center_y) + (c.z - _center_z) * (c.z - _center_z);
		double r_square = _radius * _radius;
		if (distance_square < r_square)
		{
			cell_type_sample = true;
		}

	}

	if(cell_type_sample)
	{
		enum InitShapeEnum
		{
			SHAPE_RANDOM,
			STEM_INNER,
			STEM_OUTER,
			STEM_SLICE
		};

		int initShape = params.getVal(PARAM_STEM_MODE);

		bool random = false;

		if (_stationary || initShape == InitShapeEnum::SHAPE_RANDOM) {
			random = true;
		}

		if (random) {
			int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
			type = AgentTypeEnum::CELL_TYPE_CANCER;
			state = AgentStateEnum::CANCER_PROGENITOR;

			int i = rng.sample_cdf(_celltype_cdf);
			if (i <= dmax + 1)
			{
				create_cell = true;
				if (i == 0)
				{
					state = AgentStateEnum::CANCER_STEM;
				}
				else if (i == dmax + 1)
				{
					state = AgentStateEnum::CANCER_SENESCENT;
				}
				else {
					state = AgentStateEnum::CANCER_PROGENITOR;
					div = dmax + 1 - i;
				}
			}
		}
		else { //stem arrangement cell type is determine based on its location rather than random sampled
			unsigned int dmax = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
			double csc_ratio = params.getVal(PARAM_DENSITY_CSC);
			type = AgentTypeEnum::CELL_TYPE_CANCER;
			if (initShape == InitShapeEnum::STEM_INNER) 
			{
				double r_stem = std::pow(csc_ratio, 0.5) * _radius;
				double distance_square = (c.x - _center_x) * (c.x - _center_x) + (c.y - _center_y) * (c.y - _center_y) + (c.z - _center_z) * (c.z - _center_z);
				double stem_square = r_stem * r_stem;
				if (distance_square < stem_square)
				{
					create_cell = true;
					state = AgentStateEnum::CANCER_STEM;
				}
				else {
					int i = rng.sample_cdf(_celltype_cdf);
					if (i <= dmax + 1)
					{
						create_cell = true;
						if (i == dmax + 1)
						{
							state = AgentStateEnum::CANCER_SENESCENT;
						}
						else {
							state = AgentStateEnum::CANCER_PROGENITOR;
							if (i == 0) {
								div = dmax;
							}
							else {
								div = dmax + 1 - i;
							}
						}
					}
				} 
			}
			else if (initShape == InitShapeEnum::STEM_OUTER) {
				double r_stem = std::pow(1 - csc_ratio, 0.5) * _radius;
				double distance_square = (c.x - _center_x) * (c.x - _center_x) + (c.y - _center_y) * (c.y - _center_y) + (c.z - _center_z) * (c.z - _center_z);
				double stem_square = r_stem * r_stem;
				if (distance_square > stem_square)
				{
					create_cell = true;
					state = AgentStateEnum::CANCER_STEM;
				}
				else {
					int i = rng.sample_cdf(_celltype_cdf);
					if (i <= dmax + 1)
					{
						create_cell = true;
						if (i == dmax + 1)
						{
							state = AgentStateEnum::CANCER_SENESCENT;
						}
						else {
							state = AgentStateEnum::CANCER_PROGENITOR;
							if (i == 0) {
								div = dmax;
							}
							else {
								div = dmax + 1 - i;
							}
						}
					}
				}
			}
			else if (initShape == InitShapeEnum::STEM_SLICE) {
				double pi = 3.14159265359;
				double tan_stem = std::tan(csc_ratio * 2 * pi);
				//std::cout << "tan_stem: " << tan_stem << std::endl;
				//plus 0.5 to avoid denominator == 0
				double slope = (c.z - _center_z) / (c.x - _center_x + 0.5);
				if ((slope < tan_stem) && (c.z > _center_z) && (c.x > _center_x))
				{
					create_cell = true;
					state = AgentStateEnum::CANCER_STEM;
				}
				else {
					int i = rng.sample_cdf(_celltype_cdf);
					if (i <= dmax + 1)
					{
						create_cell = true;
						if (i == dmax + 1)
						{
							state = AgentStateEnum::CANCER_SENESCENT;
						}
						else {
							state = AgentStateEnum::CANCER_PROGENITOR;
							if (i == 0) {
								div = dmax;
							}
							else {
								div = dmax + 1 - i;
							}
						}
					}
				}
			}

		}
	}
	return create_cell;
}

void VoxelContentGen::setup(bool stationary, double cancer_prob,
	int xlim, int ylim, int zlim){
	_stationary = stationary;

	_x_min= 0;
	_y_min= 0;
	_z_min= 0;
	_x_max= _x_min + xlim;
	_y_max= _y_min + ylim;
	_z_max= _z_min + zlim;

	set_cell_cdf(_stationary);
	
	return;
}

/* Fill (0,0,0), (xlim, ylim, zlim) with cancer cells
*/
void VoxelContentGen::setup(bool stationary, int xlim, int ylim, int zlim,
	int x0, int y0, int z0) {

	_stationary = stationary;

	_x_min= x0;
	_y_min= y0;
	_z_min= z0;
	_x_max= _x_min + xlim;
	_y_max= _y_min + ylim;
	_z_max= _z_min + zlim;

	set_cell_cdf(_stationary);

	/*
	double p = (1 - pstem) / double(dmax + 1);
	std::cout << "p step: " << p << std::endl;
	_celltype_cdf[0] = pstem;

	for (size_t i = 1; i < dmax+2; i++)
	{
		_celltype_cdf[i] = _celltype_cdf[i - 1] + p;
	}
	
	_celltype_cdf[dmax + 2] = 1.0;
	*/
	std::cout << "cancer cell cdf:" << std::endl;
	for (size_t i = 0; i < _celltype_cdf.size(); i++)
	{
		std::cout << i << ", " << _celltype_cdf[i] << std::endl;
	}

	return;
}

void VoxelContentGen::setup(bool stationary, double center_x, double center_y, double center_z, double radius) {

	_stationary = stationary;
	_circular = true;

	_radius = radius;
	_center_x = center_x;
	_center_y = center_y;
	_center_z = center_z;

	set_cell_cdf(_stationary);

	/*
	double p = (1 - pstem) / double(dmax + 1);
	std::cout << "p step: " << p << std::endl;
	_celltype_cdf[0] = pstem;

	for (size_t i = 1; i < dmax+2; i++)
	{
		_celltype_cdf[i] = _celltype_cdf[i - 1] + p;
	}

	_celltype_cdf[dmax + 2] = 1.0;
	*/
	std::cout << "cancer cell cdf:" << std::endl;
	for (size_t i = 0; i < _celltype_cdf.size(); i++)
	{
		std::cout << i << ", " << _celltype_cdf[i] << std::endl;
	}


	return;
}
};// end of namespace
};