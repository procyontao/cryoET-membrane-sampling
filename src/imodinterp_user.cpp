#include <algorithm>
#include <utility>
#include <vector>
#include "imodinterp.h"

void countContoursCurrObj( Iobj * const obj, int &interp, int &key, int &totalNonEmpty ) {
	totalNonEmpty = 0;
	interp = 0;
	key = 0;

	Icont *cont(nullptr);

	for(int c=0; c < csize( obj ); c++)	{
		cont = getCont( obj, c );
		if( isEmpty(cont) )
			continue;

		totalNonEmpty++;
		if( isInterpolated( cont ) )
			interp++;
		else
			key++;
	}
}


bool imod_interpolation::performInterpolation(Iobj * const out_obj, Imod * const imod) {

	Iobj * const obj  = imodObjectGet(imod);
	int interp, key, totalNonEmpty;
	countContoursCurrObj( obj, interp, key, totalNonEmpty );
	if ( key == 0 ) {
		wprint("There are no key contours to interpolated");
		return false;
	}

	int numKeyDone = 0;
	const int numContsStart = csize( obj );

	wprint("\nTCL is interpolating   0/%d @ (0%%)", key);

	for(int c=0; c < numContsStart && c < csize( obj ); c++)	{
		Icont * const cont = getCont( obj, c );
		if( isEmpty(cont) )
			continue;

		if( !isInterpolated( cont ) )
		{
			performInterpolationOnCont( obj, c, interpolation_data.interType );
			numKeyDone++;
			int percentDone = calcPercentFloor( numKeyDone, key );
			wprint("\rLYT has interpolated %d / %d @ (%d%%)", numKeyDone, key, percentDone);
		}
	}
	wprint("\r\nDone.\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n");
	{
		std::vector<std::pair<int, float>> idx_z_list;
		idx_z_list.reserve(std::size_t(obj->contsize));
		for(int no(0); no < obj->contsize; ++no) {
			idx_z_list.emplace_back(std::make_pair(no, obj->cont[no].pts[0].z));
		}
		std::sort(idx_z_list.begin(), idx_z_list.end(),
							[](const std::pair<int, float> &lhs,
							const std::pair<int, float> &rhs) -> bool {
			return lhs.second < rhs.second;
		});
		for(const auto & [idx, z] : idx_z_list) {
			Icont * const cont(imodContourDup(&obj->cont[idx]));
			imodObjectAddContour(out_obj, cont);
		}
	}

	return true;
}
