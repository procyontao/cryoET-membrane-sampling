#pragma once
#include <iosfwd>
#include <utility>
#include <imodel.h>

std::ostream & operator<<(std::ostream &os, const Ipoint &pt);

bool imod_point_inside_boundary(const Ipoint &pt, const Iobj &obj);
b3dFloat imod_point_contours_distance(const Ipoint &pt, const Iobj &obj);
Ipoint imod_point_contours_norm(const Ipoint &pt, const Iobj &obj);
Ipoint imod_point_contours_norm(const Ipoint &pt, const Icont * const cont0, const Icont * const cont1, Ipoint ref3[3]);
Imod * imod_clone_head(const Imod * const pSrc);

enum class imod_obj_op_type { sum, trail };

bool imod_obj_op(Iobj * const out, const Iobj * const lhs, const Iobj * const rhs,
								 const imod_obj_op_type op);
int imod_get_z(const Iobj * const obj, const b3dFloat z);
Icont * imod_get_zcont(const Iobj *const obj, const b3dFloat z);
Icont * imod_get_clamp_zcont(const Iobj *const obj, const b3dFloat z);
std::pair<b3dFloat, b3dFloat> imod_object_get_zbounding(const Iobj &obj);
std::pair<Ipoint, Ipoint> imod_contour_get_bounding(const Icont &cont);
std::tuple<b3dFloat, b3dFloat, b3dFloat, b3dFloat> imod_contour_get_bounding_2d(const Icont &cont);
std::pair<Ipoint, Ipoint> imod_contour_get_bounding(const Icont &lhs, const Icont &rhs);

class imod_object {
public:
    imod_object();
    imod_object(const Iobj & obj, const b3dInt32 istart, const b3dInt32 iend=0);
    imod_object(const Iobj & obj, const b3dUInt32 lkeep, const b3dUInt32 rkeep);
		imod_object(const Iobj & obj, const Iobj & ref);
    ~imod_object();
    
    imod_object(const imod_object &) = delete;
    imod_object& operator=(const imod_object&) = delete;

    operator Iobj& () const {
        return *pobj;
    }
    operator Iobj* () const {
        return pobj;
    }
protected:
    Iobj * const pobj;
};
