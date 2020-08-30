#include <cmath>
#include <cstring>

#include <algorithm>
#include <deque>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
#include <fstream>
#include <iostream>

#include <imodel.h>
#include <b3dutil.h>
#include <iimage.h>

#include <common.h>

#include "imodsup.h"
#include "imodinterp.h"

namespace __details__ {
	template<typename fT, typename iT>
	std::ptrdiff_t calculate_steps(const fT fmin, const fT fmax, const iT stepsize) {
		return std::ptrdiff_t((fmax-fmin)/stepsize + fT(.5));
	}
	template<typename fT, typename iT>
	fT grid_start_point(const fT f, const iT stepsize) {
		return std::ceil(f/stepsize)*stepsize;
	}
	namespace include {
		template<typename fT, typename iT>
		std::ptrdiff_t calculate_steps(const fT fmin, const fT fmax, const iT stepsize) {
			return ::__details__::calculate_steps(fmin, fmax, stepsize) + 2;
		}
	template<typename fT, typename iT>
		fT grid_start_point(const fT f, const iT stepsize) {
			return ::__details__::grid_start_point(f, stepsize) - stepsize;
		}		
	}
}

void imod_process(const std::string &input,
		const std::string &reference,
		const std::string &output,
		const std::size_t halfboxsize,
		const std::size_t sampling,
		const float distance_factor,
		const std::ptrdiff_t zbridge,
		const std::ptrdiff_t intmode,
		const bool debug) {
		Imod * const mAreaMod(imodRead(input.c_str()));
		if(!mAreaMod) {
			std::cout<<"Error with model file \""<<input<<"\".\n";
			return;
		}

		const Iobj * const obj(&mAreaMod->obj[0]);
		Imod * const outMod(imod_clone_head(mAreaMod));
		std::strcpy(outMod->name, "TCL's & LYT's model");
		imodNewObject(outMod);

		Imod * const refMod(imodRead(reference.c_str()));
		if(!refMod) {
			std::cout<<"Error with ref file \""<<reference<<"\".\n";
		}
		const Iobj * const refObj(&refMod->obj[0]);
		Imod * const sumMod(imodNew());
		imodNewObject(sumMod);
		Iobj * const sumObj(imodObjectGet(sumMod));
		imod_obj_op(sumObj, obj, refObj, imod_obj_op_type::sum);

		if(debug) {
			imodFileWrite(sumMod, (output+".sum.mod").c_str());
		}

		imod_interpolation ie;
		ie.interpolation_data.interType = intmodes(intmode);
		ie.interpolation_data.tilingMethod = tilingmethod::TM_AUTO;
		ie.interpolation_data.branchingMethod = branchingmethod::BR_BRANCHING_OFF;
		ie.interpolation_data.surfResolveMethod = surfacemethod::SR_CENTER_OVERLAP;
		ie.interpolation_data.ptResolveMethod = ptmethod::PT_FOUR_PTS;
		{
			//const auto [zmin, zmax] = imod_object_get_zbounding(*obj);
			//ie.interpolation_data.zBridge = int((zmax-zmin)*b3dFloat(1.25)/imodGetZScale(mAreaMod));
			ie.interpolation_data.zBridge = zbridge;
		}
		ie.interpolation_data.interSepDistBetweenConts = 160;
		ie.interpolation_data.interFractOverlap = .5f;
		ie.interpolation_data.interLineTracker = true;

		Imod * const interpSumMod(imodNew());
		imodNewObject(interpSumMod);
		Iobj * const interpSumObj(imodObjectGet(interpSumMod));
		if(debug) {
			imodFileWrite(sumMod, (output+".before.interp.mod").c_str());
		}
		ie.performInterpolation(interpSumObj, sumMod);
		if(debug) {
			imodFileWrite(interpSumMod, (output+".after.interp.mod").c_str());			
		}

		Imod * const interpMainMod(imodNew());
		imodNewObject(interpMainMod);
		Iobj * const interpMainObj(imodObjectGet(interpMainMod));
		if(debug) {
			imodFileWrite(mAreaMod, (output+".before.main.mod").c_str());			
		}
		ie.performInterpolation(interpMainObj, mAreaMod);
		if(debug) {
			imodFileWrite(interpMainMod, (output+".after.main.mod").c_str());			
		}
		std::ofstream ofs(output+".norm");

		const auto [zmin, zmax] = imod_object_get_zbounding(*interpSumObj);
		const b3dFloat fhalfboxsize(std::sqrt(3.f)*static_cast<b3dFloat>(halfboxsize)/2.f);
		const b3dFloat fDistance(distance_factor * fhalfboxsize);
		const b3dFloat boxstepsize(static_cast<b3dFloat>(2*halfboxsize)/sampling);


		const std::ptrdiff_t zsteps(__details__::include::calculate_steps(zmin, zmax, boxstepsize));
		const b3dFloat zstart(__details__::include::grid_start_point(zmin, boxstepsize));

		Ipoint ref3[3];
		Ipoint cpt;
		cpt.z = zstart;
		for(std::ptrdiff_t zstep(0); zstep < zsteps; ++zstep, cpt.z += boxstepsize) {
			Icont * const sumCont(imod_get_zcont(interpSumObj, cpt.z));
			Icont * const mainCont(imod_get_zcont(interpMainObj, cpt.z));
			Icont * const mainContA(imod_get_clamp_zcont(interpMainObj, cpt.z));

			const auto [xmin, ymin, xmax, ymax] = imod_contour_get_bounding_2d(*sumCont);
			const std::ptrdiff_t ysteps(__details__::include::calculate_steps(ymin, ymax, boxstepsize));
			const std::ptrdiff_t xsteps(__details__::include::calculate_steps(xmin, xmax, boxstepsize));
			const b3dFloat ystart(__details__::include::grid_start_point(ymin, boxstepsize));
			const b3dFloat xstart(__details__::include::grid_start_point(xmin, boxstepsize));

			bool has_points_in_current_z(false);
			cpt.y = ystart;
			for(std::ptrdiff_t ystep(0); ystep < ysteps; ++ystep, cpt.y += boxstepsize) {
				cpt.x = xstart;
				for(std::ptrdiff_t xstep(0); xstep < xsteps; ++xstep, cpt.x += boxstepsize) {
					if(imodPointInsideCont(sumCont, &cpt)) {
						continue;
					}
					int closest{0};
					const b3dFloat eDistance(imodPointContDistance(mainCont, &cpt, 1, 1, &closest));
					if(eDistance>fDistance) {
						continue;
					}
					if(!has_points_in_current_z) {
						imodNewContour(outMod);
						has_points_in_current_z = true;
					}
					imodNewPoint(outMod, &cpt);
					const Ipoint &norm(imod_point_contours_norm(cpt, mainCont, mainContA, ref3));
					ofs << cpt << " : " << norm << ", " << ref3[0] << ", " << ref3[1] << ", " << ref3[2] << "\n";
				}
			}
		}
		imodDelete(interpMainMod);
		imodDelete(interpSumMod);

		imodFileWrite(outMod, output.c_str());
		std::cout<<"Result has been written to \"" << output << "\".\n";
		imodDelete(outMod);

		imodDelete(refMod);
		imodDelete(mAreaMod);
	}

