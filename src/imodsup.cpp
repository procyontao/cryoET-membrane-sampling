#include <cstdlib>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "vec.hpp"
#include "imodsup.h"

inline Ipoint vec3_to_Ipoint(const __details__::vec3<b3dFloat> &v3) {
	return {v3.x, v3.y, v3.z};
}

inline __details__::vec3<b3dFloat> Ipoint_to_vec3(const Ipoint &pt) {
	return __details__::vec3<b3dFloat>{pt.x, pt.y, pt.z};
}

std::ostream & operator<<(std::ostream &os, const Ipoint &pt) {
	os << fmt::format("({: 9.3f},{: 9.3f},{: 9.3f})", pt.x, pt.y, pt.z) ;
	return os;
}

int imod_get_z(const Iobj * const obj, const b3dFloat z) {
	b3dFloat zDistance(std::numeric_limits<b3dFloat>::max());
	int zCo(0);
	for(int co(0); co < obj->contsize; ++co) {
		const b3dFloat cDistance(std::fabs(z-obj->cont[co].pts[0].z));
		if(zDistance>cDistance) {
			zDistance = cDistance;
			zCo = co;
		}
	}
	return zCo;
}

Icont * imod_get_zcont(const Iobj * const obj, const b3dFloat z) {
	return &obj->cont[imod_get_z(obj, z)];
}

Icont * imod_get_clamp_zcont(const Iobj *const obj, const b3dFloat z) {
	const int zCo(imod_get_z(obj, z));
	int zCoA;
	if(0 == zCo) {
		zCoA = 1;
	} else if ( obj->contsize-1 == zCo) {
		zCoA = obj->contsize-2;
	} else {
		const b3dFloat lDistance(std::fabs(z-obj->cont[zCo-1].pts[0].z));
		const b3dFloat hDistance(std::fabs(z-obj->cont[zCo+1].pts[0].z));
		if(lDistance<hDistance) {
			zCoA = zCo-1;
		} else {
			zCoA = zCo+1;
		}
	}
	return obj->cont + zCoA;
}

bool imod_obj_op(Iobj * const out, const Iobj * const lhs, const Iobj * const rhs,
								 const imod_obj_op_type op) {
	if(lhs->contsize != rhs->contsize) {
		return false;
	}
	switch(op) {
		case imod_obj_op_type::sum:
			for (int co(0); co < lhs->contsize; ++co) {
				Icont * const cont(imodContourDup(&lhs->cont[co]));
				for (int no(0); no < rhs->cont[co].psize; ++no) {
					Ipoint * const pt{&rhs->cont[co].pts[no]};
					imodPointAppend(cont, pt);
				}
				imodObjectAddContour(out, cont);
				std::free(cont);
			}
			break;
		case imod_obj_op_type::trail:
			for (int co(0); co < lhs->contsize; ++co) {
				Icont * const cont{imodContourDup(&rhs->cont[co])};
				{
					Ipoint * const pt{&lhs->cont[co].pts[lhs->cont[co].psize-1]};
					imodPointAdd(cont, pt, 0);
				}
				{
					Ipoint * const pt{&lhs->cont[co].pts[0]};
					imodPointAdd(cont, pt, 0);
				}
				imodObjectAddContour(out, cont);
				std::free(cont);
			}
			break;
	};
	return true;
}


bool imod_point_inside_boundary(const Ipoint &pt, const Iobj &obj) {
	b3dFloat min_z{std::numeric_limits<b3dFloat>::max()};
	b3dFloat max_z{-std::numeric_limits<b3dFloat>::max()};
	for (int co(0); co < obj.contsize; ++co) {
		if (obj.cont[co].psize) {
			const b3dFloat cur_z{obj.cont[co].pts[0].z};
			if(min_z > cur_z) {
				min_z = cur_z;
			}
			if(max_z < cur_z) {
				max_z = cur_z;
			}
		}
	}
	if(pt.z < min_z || pt.z > max_z) {
		return false;
	}
	b3dFloat min_distance{std::numeric_limits<b3dFloat>::max()};
	int no{-1};
	for (int co(0); co < obj.contsize; ++co) {
		if (obj.cont[co].psize) {
			const b3dFloat cur_distance{std::abs(obj.cont[co].pts[0].z - pt.z)};
			if(cur_distance < min_distance) {
				min_distance = cur_distance;
				no = co;
			}
		}
	}
	if(imodPointInsideCont(&obj.cont[no], const_cast<Ipoint *>(&pt))) {
		return true;
	} else {
		return false;
	}
}

inline b3dFloat imod_point_contour_2nd_distance(Ipoint * const pt, Icont * const pCont, const int closest, Ipoint &pOut) {
	b3dFloat distance(std::numeric_limits<b3dFloat>::max());
	int no_2nd;
	for(int no(0); no < pCont->psize && no != closest; ++no) {
		const b3dFloat cur_distance(imodPointDistance(pt, &pCont->pts[no]));
		if(cur_distance < distance) {
			distance = cur_distance;
			no_2nd = no;
		}
	}
	pOut = pCont->pts[no_2nd];
	return distance;
}

Ipoint imod_point_contours_norm(const Ipoint &pt, const Iobj &obj) {
	b3dFloat distance_0(std::numeric_limits<b3dFloat>::max());
	int co_id_0, closest_0;
	int closest;
	for (int co(0); co < obj.contsize; ++co) {
		const b3dFloat d(imodPointContDistance(&obj.cont[co], const_cast<Ipoint*const>(&pt), 1, 1, &closest));
		if(d < distance_0) { distance_0 = d; co_id_0 = co; closest_0 = closest;}
	}
	const __details__::vec3<b3dFloat> pt0(Ipoint_to_vec3(obj.cont[co_id_0].pts[closest_0]));

	b3dFloat distance_1;
	int co_id_1, closest_1;
	if(0 == co_id_0 || obj.contsize-1 == co_id_0) {
		if(0 == co_id_0) {
			co_id_1 = co_id_0 + 1;
		} else {
			co_id_1 = co_id_0 - 1;
		}
		distance_1 = imodPointContDistance(&obj.cont[co_id_1], const_cast<Ipoint*const>(&pt), 1, 1, &closest_1);
	} else {
		int closest_l, closest_h;
		const b3dFloat distance_l(imodPointContDistance(&obj.cont[co_id_0-1], const_cast<Ipoint*const>(&pt), 1, 1, &closest_l));
		const b3dFloat distance_h(imodPointContDistance(&obj.cont[co_id_0+1], const_cast<Ipoint*const>(&pt), 1, 1, &closest_h));
		if(distance_l < distance_h) {
			distance_1 = distance_l;
			co_id_1 = co_id_0 - 1;
			closest_1 = closest_l;
		} else {
			distance_1 = distance_h;
			co_id_1 = co_id_0 + 1;
			closest_1 = closest_h;
		}
	}
	const __details__::vec3<b3dFloat> pt1(Ipoint_to_vec3(obj.cont[co_id_1].pts[closest_1]));
	__details__::vec3<b3dFloat> pt2; {
		Ipoint point0a, point1a;
		const b3dFloat distance_0a(imod_point_contour_2nd_distance(const_cast<Ipoint*const>(&pt), &obj.cont[co_id_0], closest_0, point0a));
		const b3dFloat distance_1a(imod_point_contour_2nd_distance(const_cast<Ipoint*const>(&pt), &obj.cont[co_id_1], closest_1, point1a));
		if(distance_0a < distance_1a) {
			pt2 = Ipoint_to_vec3(point0a);
		} else {
			pt2 = Ipoint_to_vec3(point1a);
		}
	}
	const __details__::vec3<b3dFloat> ptA(pt1-pt0);
	const __details__::vec3<b3dFloat> ptB(pt2-pt0);
	__details__::vec3<b3dFloat> norm(__details__::cross(ptA, ptB));
	const __details__::vec3<b3dFloat> ptr(Ipoint_to_vec3(pt));
	if(__details__::dot(norm, ptr)<b3dFloat(0)) {
		norm *= b3dFloat(-1);
	}
	return vec3_to_Ipoint(__details__::normalize(norm));
}

inline int imod_point_contours_2nd_closest(const Icont * const cont, const int closest, const Ipoint &pt) {
	int closetA;
	if(0 == closest) {
		closetA = 1;
	} else if (cont->psize-1 == closest) {
		closetA = cont->psize-2;
	} else {
		const b3dFloat lDistance(imodPointDistance(const_cast<Ipoint*const>(&pt), cont->pts+closest-1));
		const b3dFloat rDistance(imodPointDistance(const_cast<Ipoint*const>(&pt), cont->pts+closest+1));
		if(lDistance<rDistance) {
			closetA = closest - 1;
		} else {
			closetA = closest + 1;
		}
	}
	return closetA;
}

Ipoint imod_point_contours_norm(const Ipoint &pt, const Icont * const cont0, const Icont * const cont1, Ipoint ref3[3]) {
	int closest_l, closest_h;
	imodPointContDistance(const_cast<Icont*const>(cont0), const_cast<Ipoint*const>(&pt), 1, 1, &closest_l);
	imodPointContDistance(const_cast<Icont*const>(cont1), const_cast<Ipoint*const>(&pt), 1, 1, &closest_h);
	const Ipoint * ppt3(nullptr); {
		const Ipoint &ptL(cont0->pts[imod_point_contours_2nd_closest(const_cast<Icont*const>(cont0), closest_l, pt)]);
		const Ipoint &ptH(cont1->pts[imod_point_contours_2nd_closest(const_cast<Icont*const>(cont1), closest_h, pt)]);
		const b3dFloat dL(imodPointDistance(const_cast<Ipoint*const>(&pt), const_cast<Ipoint*const>(&ptL)));
		const b3dFloat dH(imodPointDistance(const_cast<Ipoint*const>(&pt), const_cast<Ipoint*const>(&ptH)));
		if(dL < dH) {
			ppt3 = &ptL;
		} else {
			ppt3 = &ptH;
		}
	}

	ref3[0] = cont0->pts[closest_l];
	ref3[1] = cont1->pts[closest_h];
	ref3[2] = *ppt3;

	const __details__::vec3<b3dFloat> pt0(Ipoint_to_vec3(cont0->pts[closest_l]));
	const __details__::vec3<b3dFloat> pt1(Ipoint_to_vec3(cont1->pts[closest_h]));
	const __details__::vec3<b3dFloat> pt2(Ipoint_to_vec3(*ppt3));

	const __details__::vec3<b3dFloat> ptA(pt1-pt0);
	const __details__::vec3<b3dFloat> ptB(pt2-pt0);
	__details__::vec3<b3dFloat> norm(__details__::cross(ptA, ptB));
	const __details__::vec3<b3dFloat> ptr(Ipoint_to_vec3(pt)-pt0);
	if(__details__::dot(norm, ptr)<b3dFloat(0)) {
		norm *= b3dFloat(-1);
	}
	return vec3_to_Ipoint(__details__::normalize(norm));
}

b3dFloat imod_point_contours_distance(const Ipoint &pt, const Iobj &obj) {
	b3dFloat distance(std::numeric_limits<b3dFloat>::max());
	int closest;
	for (int co(0); co < obj.contsize; ++co) {
		const b3dFloat d(imodPointContDistance(&obj.cont[co], const_cast<Ipoint*const>(&pt), 1, 1, &closest));
		if(d < distance) { distance = d; }
	}
	return distance;
}

std::pair<b3dFloat, b3dFloat> imod_object_get_zbounding(const Iobj &obj) {
	b3dFloat zmin(std::numeric_limits<b3dFloat>::max());
	b3dFloat zmax(-std::numeric_limits<b3dFloat>::max());
	for(int co(0); co < obj.contsize; ++co) {
		for(int no(0); no < obj.cont[co].psize; ++no) {
			const b3dFloat cz(obj.cont[co].pts[no].z);
			if(zmin>cz) { zmin = cz; }
			if(zmax<cz) { zmax = cz; }
		}
	}
	return std::make_pair(zmin, zmax);
}

std::pair<Ipoint, Ipoint> imod_contour_get_bounding(const Icont &cont) {
	Ipoint pmin{std::numeric_limits<b3dFloat>::max(), std::numeric_limits<b3dFloat>::max(), std::numeric_limits<b3dFloat>::max()};
	Ipoint pmax{-std::numeric_limits<b3dFloat>::max(), -std::numeric_limits<b3dFloat>::max(),-std::numeric_limits<b3dFloat>::max()};
	for(int no(0); no < cont.psize; ++no) {
		const Ipoint & pt(cont.pts[no]);
		if(pmin.x>pt.x) { pmin.x = pt.x; }
		if(pmin.y>pt.y) { pmin.y = pt.y; }
		if(pmin.z>pt.z) { pmin.z = pt.z; }
		if(pmax.x<pt.x) { pmax.x = pt.x; }
		if(pmax.y<pt.y) { pmax.y = pt.y; }
		if(pmax.z<pt.z) { pmax.z = pt.z; }
	}
	return std::make_pair(pmin, pmax);
}

std::tuple<b3dFloat, b3dFloat, b3dFloat, b3dFloat> imod_contour_get_bounding_2d(const Icont &cont) {
	b3dFloat xmin{std::numeric_limits<b3dFloat>::max()}, ymin{std::numeric_limits<b3dFloat>::max()};
	b3dFloat xmax{-std::numeric_limits<b3dFloat>::max()}, ymax{-std::numeric_limits<b3dFloat>::max()};
	for(int no(0); no < cont.psize; ++no) {
		const Ipoint & pt(cont.pts[no]);
		if(xmin>pt.x) { xmin = pt.x; }
		if(ymin>pt.y) { ymin = pt.y; }
		if(xmax<pt.x) { xmax = pt.x; }
		if(ymax<pt.y) { ymax = pt.y; }
	}
	return std::make_tuple(xmin, ymin, xmax, ymax);
}

std::pair<Ipoint, Ipoint> imod_contour_get_bounding(const Icont &lhs, const Icont &rhs) {
	auto [lmin, lmax] = imod_contour_get_bounding(lhs);
	auto [rmin, rmax] = imod_contour_get_bounding(rhs);
	return std::make_pair<Ipoint, Ipoint>(
	{std::min(lmin.x, rmin.x), std::min(lmin.y, rmin.y), std::min(lmin.z, rmin.z)},
	{std::max(lmax.x, rmax.x), std::max(lmax.x, rmax.x), std::max(lmax.x, rmax.z)});
}


Imod * imod_clone_head(const Imod * const pSrc) {
	Imod * const pDest(imodNew());
	pDest->flags = pSrc->flags;
	pDest->drawmode = pSrc->drawmode;
	pDest->mousemode = pSrc->mousemode;
	pDest->blacklevel = pSrc->blacklevel;
	pDest->whitelevel = pSrc->whitelevel;
	pDest->xoffset = pSrc->xoffset;
	pDest->yoffset = pSrc->yoffset;
	pDest->zoffset = pSrc->zoffset;
	pDest->xscale = pSrc->xscale;
	pDest->yscale = pSrc->yscale;
	pDest->zscale = pSrc->zscale;
	pDest->pixsize = pSrc->pixsize;
	pDest->units = pSrc->units;
	pDest->alpha = pSrc->alpha;
	pDest->beta = pSrc->beta;
	pDest->gamma = pSrc->gamma;
	pDest->xybin = pSrc->xybin;
	pDest->zbin = pSrc->zbin;
	return pDest;
}

imod_object::imod_object():pobj(imodObjectNew()) {}

imod_object::imod_object(const Iobj & obj, const b3dInt32 istart, const b3dInt32 iend) 
    :pobj(imodObjectNew()) {
    for (b3dInt32 co(0); co < obj.contsize; ++co) {
        const b3dInt32 psize(obj.cont[co].psize);
		if ( psize < 1 ||
             std::abs(istart) >= psize ||
             std::abs(iend) > psize ) {
            continue;
        }
        const b3dInt32 pstart(istart<0?istart+psize:istart);
        const b3dInt32 pend(iend<=0?iend+psize:iend);
        if(pstart>=pend) {
            continue;
        }
        Icont * const pcont(imodContourNew());
        for(b3dInt32 no(pstart); no<pend; ++no) {
            imodPointAppend(pcont, &obj.cont[co].pts[no]);
        }
        imodObjectAddContour(pobj, pcont);
        std::free(pcont);
    }
}

imod_object::imod_object(const Iobj & obj, const b3dUInt32 lkeep, const b3dUInt32 rkeep) 
    :pobj(imodObjectNew()) {
    if(lkeep+rkeep>0) {
        for (b3dInt32 co(0); co < obj.contsize; ++co) {
						const b3dUInt32 psize(static_cast<b3dUInt32>(obj.cont[co].psize));
            if ( psize < 1 ||
                lkeep + rkeep > psize) {
                continue;
            }
            Icont * const pcont(imodContourNew());
			for(b3dUInt32 no(0); no<lkeep; ++no) {
                imodPointAppend(pcont, &obj.cont[co].pts[no]);
            }
            for(b3dUInt32 no(psize-rkeep); no<psize; ++no) {
                imodPointAppend(pcont, &obj.cont[co].pts[no]);
            }
            imodObjectAddContour(pobj, pcont);
            std::free(pcont);
        }
    }
}

/*
imod_object::imod_object(const Iobj & obj, const Iobj & ref)
	:pobj(imodObjectNew()) {

}
*/

imod_object::~imod_object() {
	imodObjectDelete(pobj);
}
