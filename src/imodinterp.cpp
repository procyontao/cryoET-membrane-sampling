#include <cmath>
#include <climits>
#include <deque>
#include <limits>

#include <imodel.h>
#include <icont.h>
#include <ipoint.h>
#include <imesh.h>
#include <mkmesh.h>

#include "imodsup.h"
#include "an_common.h"
#include "an_cont.h"
#include "imodinterp.h"

void ContInfo::resetAll()	{
		idx = -1;
		setPt( &ll, -1, -1, -1 );
		setPt( &ur, -1, -1, -1 );
		setPt( &centerPt, -1, -1, -1 );
		z = -1;
		loaded  = false;
		checked = 0;
}

//-------------------------------
//** CONSTRUCTORS:



imod_interpolation::imod_interpolation()
	:obj(nullptr) {

}


//------------------------
//-- Resets all variables in the imod_interpolation structure to defaults,
//-- clears the contour information vector and resizes/empies both z tables.

void imod_interpolation::resetAll() {
	obj      = nullptr;
	sContIdx = -1;
	sContZ   = -1;

	conti.clear();

	if( !ztableKey.empty() )
		for(std::size_t i=0; i<ztableKey.size(); i++ )
			ztableKey[i].idxs.clear();

	if( !ztableInt.empty() )
		for(std::size_t i=0; i<ztableInt.size(); i++ )
			ztableInt[i].idxs.clear();

	closestAboveIdx = -1;
	closestBelowIdx = -1;

	aboveIdxs.clear();
	belowIdxs.clear();
}



//------------------------
//-- Is typically executed after a new (key) contour is added/modified by
//-- the user and we want to perform interpolation on it. This function
//-- first tries to find the nearest key contour above and below the
//-- new one (within a set range), and, depending on if any are found
//-- (and also the type of interpolation selected) the function will
//-- delete nearby interpolated contours in the way (above and below) and
//-- generate series of new interpolated contours.

void imod_interpolation::performInterpolationOnCont( Iobj *_obj, int _contIdx, const intmodes interpolationType,
																										 const float modScaleZ) {
	//## SETUP INTERPOLATION EVENT DATA (POINTERS, VECTORS, Z TABLES, ETC):
	resetAll();

	obj     = _obj;
	sContIdx = _contIdx;
	sContZ   = getZInt( getC(sContIdx) );

	regenerateContInfoVector();
	regenerateZTables();
	regenerateConnections();

	//## PERFORM INTERPOLATION:

	int numAdded = 0;
	int numFlagged = 0;
	int numDeleted = 0;

	numFlagged = deleteInterpolatedContsBetweenKeyContsAboveAndBelow();

	switch( interpolationType )	{
		case (intmodes::INT_NO_INTERPOLATE):
			break;

		case (intmodes::INT_LINEAR):
			numAdded = interp_Linear( sContIdx, interpolation_data.zBridge );
			break;

		case (intmodes::INT_SPHERICAL):
			interp_Spherical( sContIdx, interpolation_data.zBridge, modScaleZ );
			break;

		case (intmodes::INT_SMOOTH):
			interp_SmoothCrudeAll( sContIdx, interpolation_data.zBridge );
			break;

		case (intmodes::INT_SMOOTH_PTS):
			interp_SmoothPointwise( sContIdx, interpolation_data.zBridge );
			break;

		default:
			wprint("This type of interpolation not implemented yet\n");
			break;
	}

	numDeleted = removeAllDeleteFlaggedContoursFromObj( obj );

	if(numDeleted < numFlagged)
		wprint( "ERROR - not all conts flagged for delete were deleted\n" );

	//wprint("interpolation cont %d, deleted %d, added %d\n",sContIdx,numDeleted,numAdded);
}


//------------------------
//-- Find the nearest key contour above and below the contour given,
//-- then deletes all the interpolated contours in between which appear to
//-- belong to the same surface

int imod_interpolation::
deleteImmediatelyAdjacentInterpolatedConts( Iobj *_obj, int _contIdx )
{
	//## SETUP INTERPOLATION EVENT DATA (POINTERS, VECTORS, Z TABLES, ETC):
	resetAll();

	obj     = _obj;
	sContIdx = _contIdx;
	sContZ   = getZInt( getC(sContIdx) );

	regenerateContInfoVector();
	regenerateZTables();
	regenerateConnections();

	//## PERFORM DELETION:

	deleteInterpolatedContsBetweenKeyContsAboveAndBelow();
	return ( removeAllDeleteFlaggedContoursFromObj( obj ) );
}


//------------------------
//-- Attempts to find a list of consequtive (same-surface) interpolated contours
//-- with greater than "minSpanSize" contours, then returns the index of the
//-- one in the middle. If none is found it returns -1.
//-- The search begins from the contour in "_obj" at "_contIdx".

int imod_interpolation::
findMiddleNextLargeInterpolatedSpan( Iobj *_obj, int _contIdx, int minSpanSize )
{
	//## SETUP INTERPOLATION EVENT DATA (POINTERS, VECTORS, Z TABLES, ETC):
	resetAll();

	obj     = _obj;
	sContIdx = _contIdx;
	sContZ   = getZInt( getC(sContIdx) );

	regenerateContInfoVector();
	regenerateZTables();

	//## PERFORM DELETION:

	for(int i=0; i<int(conti.size()); i++ )     // for each interpolated contour (contI):
	{
		int cIdxI = (i+sContIdx+1) % csize(obj);
		Icont *contI = getC(cIdxI);
		if( !isInterpolated(contI) || conti[std::size_t(cIdxI)].checked || psize(contI)==0 )
			continue;

		if( cIdxI == 0 )
			wprint("Starting from beginning ...\n");

		int contZ = getCZ(cIdxI);

		//## COMPILE LIST OF CONNECTED INTERPOLATED CONTOURS ABOVE AND BELOW CONTI:

		std::vector<int> interpConts;            // list of connected, interpolated contours
		interpConts.resize(0);
		interpConts.push_back(cIdxI);
		conti[cIdxI].checked = 1;
		int j, cIdxPrev;

		for(int z=contZ+1; z<interpolation_data.zsize; z++)  // for each slice above current contour:
		{
			for(j=0; j<numIntContsAtZ(z); j++)  // for each key contour at this z:
			{
				int cIdxJ = idxIntContZ(z,j);
				if ( conti[cIdxJ].checked == 0
						 && contoursSameSurf(cIdxJ,cIdxI) )	{ // if belongs to same surface:
					interpConts.push_back( cIdxJ );
					conti[cIdxJ].checked = 1;
					break;
				}
			}
			if(j==numIntContsAtZ(z))     // if no connected cont found at this z: break
				break;
		}
		for(int z=contZ-1; z>=0; z--)           // for each slice below current contour:
		{
			for(j=0; j<numIntContsAtZ(z); j++)  // for each key contour at this z:
			{
				int cIdxJ = idxIntContZ(z,j);
				if ( conti[cIdxJ].checked == 0
						 && contoursSameSurf(cIdxJ,cIdxI) )	{ // if belongs to same surface:
					interpConts.insert( interpConts.begin(), cIdxJ );
					conti[cIdxJ].checked = 1;
					break;
				}
			}
			if(i==numIntContsAtZ(z))     // if no connected cont found at this z: break
				break;
		}

		int numInterpConts = (int)interpConts.size();    // number of connected/consequtive
																										 // interpoalted contours found

		if( numInterpConts >= minSpanSize )
		{
			int middleMostContIdx = interpConts[ int(numInterpConts/2) ];
			wprint("Span of %d interpolated contours found\n", numInterpConts);
			return ( middleMostContIdx );
		}
	}

	return -1;
}


//------------------------
//-- Resizes and populates the contour information vector (conti),
//-- calculating the minimum bounding rectangle and other information
//-- for any contours which lie in the given Z range


void imod_interpolation::
regenerateContInfoVector( int _minZLimit, int _maxZLimit, int startIdx )
{
	minZLimit = _minZLimit;
	maxZLimit = _maxZLimit;

	//## POPULATE CONTOURINFO VECTOR:

	conti.resize( csize(obj) );
	for(int c=startIdx; c<csize(obj); c++ )
	{
		//conti[c].resetAll();

		Icont *cont = getC(c);
		conti[c].idx = c;
		conti[c].z   = getZInt(cont);

		//if( conti[c].z < minZLimit || conti[c].z > maxZLimit )
		//  continue;

		cont_getMBR( cont, &conti[c].ll, &conti[c].ur );
		conti[c].centerPt = line_getPtHalfwayBetween( &conti[c].ll, &conti[c].ur );
		conti[c].loaded = true;
	}
}


//------------------------
//-- Regenerates the Z table for key contours and interpolated contours,
//-- which list the index of contours (within the object) on each z slice.

void imod_interpolation::regenerateZTables(  )
{
	//## POPULATE Z TABLES:

	ztableKey.resize( interpolation_data.zsize );
	ztableInt.resize( interpolation_data.zsize );

	for(int c=0; c<csize(obj); c++ )
	{
		Icont *cont = getC(c);

		if( isEmpty(cont) )
			continue;

		int z = getZInt(cont);
		if( z<0 || z>=interpolation_data.zsize )
			continue;

		if( isInterpolated(cont) )
			ztableInt[z].idxs.push_back(c);
		else
			ztableKey[z].idxs.push_back(c);
	}
}

//------------------------
//-- Redetermines the closest contour(s) above and below the contour
//-- currently being interpolated.
//-- Input:  'sContIdx'
//-- Output: 'aboveIdxs', 'belowIdxs', 'closestAboveIdx', 'closestBelowIdx'

void imod_interpolation::regenerateConnections(  )
{
	//## FIND CLOSEST KEY CONTOURS (CLOSEST IN Z) WHICH SEEM TO BELONG TO THE
	//## SAME SURFACE ABOVE AND BELOW THE BASE CONTOUR:

	aboveIdxs = findIdxsNearestKeyConts(sContIdx,true,interpolation_data.zBridge);
	belowIdxs = findIdxsNearestKeyConts(sContIdx,false,interpolation_data.zBridge);

	//## FIND CLOSEST KEY CONTOURS IN X & Y FROM THE LISTS ABOVE:

	closestAboveIdx = findIdxClosestKeyContInList(sContIdx,aboveIdxs);
	closestBelowIdx = findIdxClosestKeyContInList(sContIdx,belowIdxs);
}


//------------------------
//-- Flags interpolated contours which are between the start contour (sContIdx)
//-- and the closest key contour above and below (and appear to be in the same
//-- surface) for deletion.

int imod_interpolation::deleteInterpolatedContsBetweenKeyContsAboveAndBelow()
{
		int deleteRangeMaxZ = isKeyContAbove() ? getCZ(closestAboveIdx) : sContZ;
		int deleteRangeMinZ = isKeyContBelow() ? getCZ(closestBelowIdx) : sContZ;
		return deleteAllSameSurfInterpContsInZRange(deleteRangeMinZ,deleteRangeMaxZ,sContIdx);
}


//-------------------------------
//** SURFACE RESOLVING METHODS:





//------------------------
//-- Determines if two contours (appear to) belong to the same surface.
//-- It (typically) does this by determining the amount by which they overlap.

bool imod_interpolation::
contoursSameSurf( int c1Idx, int c2Idx )
{
	Icont *cont1 = getC(c1Idx);
	Icont *cont2 = getC(c2Idx);

	const bool mbrsTouch = mbr_doBBoxesOverlap2D( getLL(c1Idx), getUR(c1Idx),
																					getLL(c2Idx), getUR(c2Idx) );

	switch (interpolation_data.surfResolveMethod)
	{
		case( surfacemethod::SR_TOUCHING ):
		{
			return ( mbrsTouch && cont_doContsTouch(cont1,cont2) );
		}
			break;

		case( surfacemethod::SR_CENTER_OVERLAP ):
		{
			if( !mbrsTouch )
				return false;
			return (   imodPointInsideCont(cont1,getCenterPt(c2Idx))
							|| imodPointInsideCont(cont2,getCenterPt(c1Idx)) );
		}
			break;

		case( surfacemethod::SR_WITHIN_MIN_DIST ):
		{
			float distBetweenMBRs = mbr_distBetweenBBoxes2D( getLL(c1Idx), getUR(c1Idx),
																											 getLL(c2Idx), getUR(c2Idx) );
			if( distBetweenMBRs > interpolation_data.interSepDistBetweenConts )
				return false;
			float minDist = cont_minDistBetweenContPts2D(cont1, cont2, true);
			return ( minDist <= interpolation_data.interSepDistBetweenConts );
		}
			break;

		case( surfacemethod::SR_FRACTION_OVERLAP ):
		{
			if( !mbrsTouch )
				return false;

			float frac1, fract2;
			Icont *cs1 = imodel_contour_scan( cont1 );
			Icont *cs2 = imodel_contour_scan( cont2 );
			imodel_overlap_fractions( &cs1,*getLL(c1Idx),*getUR(c1Idx),
																&cs2,*getLL(c2Idx),*getUR(c2Idx),
																&frac1,&fract2 );
			imodContourDelete( cs1 );
			imodContourDelete( cs2 );

			float fractionOverlap = std::max( frac1, fract2 );
			return ( fractionOverlap >= interpolation_data.interFractOverlap );
		}
			break;

		case( surfacemethod::SR_MBR_TOUCH ):
		{
			return (mbrsTouch);
		}
			break;



		default:
			return (false);
	}
}


//------------------------
//-- Determines if two contours (appear to) belong to the same surface.
//-- It (typically) does this by determining the amount by which they overlap.

bool imod_interpolation::
contoursSameSurfSlow( Icont *cont1, Icont *cont2 )
{
	Ipoint cont1ll, cont1ur;
	Ipoint cont2ll, cont2ur;

	imodContourGetBBox( cont1, &cont1ll, &cont1ur );
	imodContourGetBBox( cont2, &cont2ll, &cont2ur );

	Ipoint centerMBRCont1 = line_getPtHalfwayBetween( &cont1ll, &cont1ur );
	Ipoint centerMBRCont2 = line_getPtHalfwayBetween( &cont2ll, &cont2ur );

	const bool mbrsTouch = mbr_doBBoxesOverlap2D( &cont1ll, &cont1ur,
																					&cont2ll, &cont2ur );


	switch (interpolation_data.surfResolveMethod)
	{
		case( surfacemethod::SR_TOUCHING ):
		{
			return ( mbrsTouch && cont_doContsTouch(cont1,cont2) );
		}
			break;

		case( surfacemethod::SR_CENTER_OVERLAP ):
		{
			if( !mbrsTouch )
				return false;
			return (   imodPointInsideCont(cont1,&centerMBRCont2)
								 || imodPointInsideCont(cont2,&centerMBRCont1) );
		}
			break;

		case( surfacemethod::SR_WITHIN_MIN_DIST ):
		{
			float distBetweenMBRs = mbr_distBetweenBBoxes2D( &cont1ll, &cont1ur,
																											 &cont2ll, &cont2ur );
			if( distBetweenMBRs > interpolation_data.interSepDistBetweenConts )
				return false;
			float minDist = cont_minDistBetweenContPts2D(cont1, cont2, true);
			return ( minDist <= interpolation_data.interSepDistBetweenConts );
		}
			break;

		case( surfacemethod::SR_FRACTION_OVERLAP ):
		{
			if( !mbrsTouch )
				return false;

			float frac1, fract2;
			Icont *cs1 = imodel_contour_scan( cont1 );
			Icont *cs2 = imodel_contour_scan( cont2 );
			imodel_overlap_fractions( &cs1,cont1ll,cont1ur,
																&cs2,cont2ll,cont2ur,
																&frac1,&fract2 );
			imodContourDelete( cs1 );
			imodContourDelete( cs2 );

			float fractionOverlap = std::max( frac1, fract2 );
			return ( fractionOverlap >= interpolation_data.interFractOverlap );
		}
			break;

		case( surfacemethod::SR_MBR_TOUCH ):
		{
			return (mbrsTouch);
		}
			break;



		default:
			return (false);
	}
}

//------------------------
//-- Takes a base contour and finds/returns the indexes of the nearest key
//-- contours (nearest in z) which (seem to) belong to the same surface and
//-- are within a specified range of z slices relative to the base contours).
//-- NOTE: Even if "minZDist" is 0, surfaces on the SAME z slice are skipped.
//-- NOTE: By changing "maxDistAway" to "minZRelativeSlice" and using a
//--       "maxZVal" and "zDistFromBase" integers can easily change this function
//--       to search both ABOVE and BELOW the base contour - but haven't found
//--       a need yet.

std::vector<int> imod_interpolation::
findIdxsNearestKeyConts( int cbIdx, bool above, int maxZDist, int minZDist  )
{
	if( minZDist < 0 || maxZDist < 0 )	{   // should never happen
		wprint( "ERROR/MISUSE: findIdxsNearestKeyConts()" );
	}

	std::vector<int> closestContIdxs;		// stores indexes of closest (same-surface) conts found
	int baseContZ = getCZ( cbIdx );
	int lastZLists = (int)ztableKey.size() - 1;

	if(above)
	{
		int minZ  = MAX( baseContZ + minZDist, 0 );
		int maxZ  = MIN( baseContZ + maxZDist, lastZLists );

		for(int z=minZ; z<=maxZ; z++)        // for each z list in range:
		{
			for(int i=0; i<numKeyContsAtZ(z); i++)  // for each key contour at this z:
			{
				int cidx = idxKeyContZ(z,i);
				if ( contoursSameSurf(cidx,cbIdx) )	 // if belongs to same surface:
					closestContIdxs.push_back( cidx );			// add it to list of "closest conts"
			}
			if( closestContIdxs.size() > 0 )    // if matching contours were found:
				break;                              // return them!
		}
	}
	else
	{
		int minZ  = MIN( baseContZ - minZDist, lastZLists );
		int maxZ  = MAX( baseContZ - maxZDist, 0 );

		for(int z=minZ; z>=maxZ; z--)       // for each z list in range:
		{
			for(int i=0; i<numKeyContsAtZ(z); i++)  // for each key contour at this z:
			{
				int cidx = idxKeyContZ(z,i);
				if ( contoursSameSurf(cidx,cbIdx) )	 // if belongs to same surface:
					closestContIdxs.push_back( cidx );			// add it to list of "closest conts"
			}
			if( closestContIdxs.size() > 0 )    // if matching contours were found:
				break;                              // return them!
		}
	}

	return closestContIdxs;
}



//------------------------
//-- Takes a base contour and finds the nearest key contour based on centroid.
//-- Returns -1 if vector is empty.

int imod_interpolation::
findIdxClosestKeyContInList( int cbIdx, std::vector<int> conts )
{
	if(conts.empty())             // if no contours were found: return -1
	{
		return -1;
	}
	else if (conts.size() == 1)   // if one closest contour was found: return it's index
	{
		return conts[0];
	}
	else                    // else (if multiple contours were found on the same z slice):
	{                       // find & return the one has the closest centroid.

		float closestDistToCentroid = FLOAT_MAX;
		int   closestContIdx = 0;

		for(int i=0; i<int(conts.size()); i++)
		{
			float distToThisCentroid = imodel_point_dist( getCenterPt( conts[i] ),
																										getCenterPt( cbIdx) );
			if ( distToThisCentroid < closestDistToCentroid )
			{
				closestDistToCentroid = distToThisCentroid;
				closestContIdx = conts[i];
			}
		}
		return closestContIdx;
	}
}


//------------------------
//-- Takes a base contour and finds the nearest key contour above or below it, within
//-- a specified number of slices, which (seems to) belong to the same surface.
//-- If no contour matches this descritpion -1 is returned.
//-- If multiple overlapping contours on the same Z slices are found, the one with
//-- the CLOSEST centroid to the base contour is returned.

int imod_interpolation::
findIdxNearestKeyCont( int cbIdx, int maxDist, bool above )
{
	std::vector<int> closestConts = findIdxsNearestKeyConts(cbIdx,above,maxDist);
	return findIdxClosestKeyContInList( cbIdx, closestConts );
}

//------------------------
//-- Takes a base contour and finds/returns the indexes of all key
//-- contours which appear to belong to the same surface, by searching
//-- iteratively down and then up. This list will be ordered from lowest
//-- to highest contour (base contour included) HOWEVER ignores any branching:
//-- if the contour splits into two contours, only the nearest
//-- will be considered.

std::vector<int> imod_interpolation::
findIdxsAllKeyContsNoBranching( int idxBaseCont, int maxDist )
{
	std::vector<int> contsInSurface;

	//## SEARCH DOWN:
	int idxCurrCont = idxBaseCont;
	while (true) {
		int idxNextCont = findIdxNearestKeyCont( idxCurrCont, maxDist, false );
																								// see if there is a contour below this
		if ( idxNextCont == -1 )                    // if no contour if found: break
			break;

		contsInSurface.insert( contsInSurface.begin(), idxNextCont );
																			// add newly found contour to start of list
		idxCurrCont = idxNextCont;        // make this the current contour - (will
																			// check if one above next iteration)
	}

	contsInSurface.push_back( idxBaseCont );     // add base contour's idx (in the middle)


	//## SEARCH UP:
	idxCurrCont = idxBaseCont;
	while (true) {
		int idxNextCont = findIdxNearestKeyCont( idxCurrCont, maxDist, true );
																							// see if there is a contour above this
		if ( idxNextCont == -1 )									// if no contour if found: break
			break;

		contsInSurface.push_back(idxNextCont);		// add newly found contour to end of list
		idxCurrCont = idxNextCont;								// make this the current contour - (will
																							// check if one above next iteration)
	}

	return contsInSurface;
}

//------------------------
//-- Deletes any interpolated contours (in the model) which are in the given range
//-- of slices and appears to belong to the same surface as the contour provided.

int imod_interpolation::
deleteAllSameSurfInterpContsInZRange( int minZ, int maxZ, int cbIdx )
{
	int numDeleted = 0;

	minZ = MAX( minZ, 0 );
	maxZ = MIN( maxZ,interpolation_data.zsize-1 );

	for(int z=minZ; z<=maxZ; z++)
	{
		for(int i=0; i<numIntContsAtZ(z); i++)
		{
			int cidx = idxIntContZ(z,i);
			if ( contoursSameSurf(cidx,cbIdx) )	 // if belongs to same surface:
			{
				setDeleteFlag( getC(cidx), 1 );
				numDeleted++;
			}
		}
	}
	return numDeleted;
}

//------------------------
//-- Deletes all adjacent interpolated contours (in the model) above and below
//-- the given contour. Unlike "deleteAllSameSurfInterpContsInZRange",
//-- this works by checking each adjacent slice in order, finding a matching
//-- interpolated contour, then using this to check the next slice until no more
//-- matching interpolated contour is found... (no range is specified).

int imod_interpolation::
deleteAllSameSurfInterpContsEitherSideOfCont( int cbIdx )
{
	int numDeleted = 0;

	//## SEARCH DOWN AND FLAG MATCHING INTERPOLATED CONTOURS FOR DELETION:

	int prevFoundIdx = cbIdx;              // countour to check for overlap
	bool matchFound = true;
	for (int z=getCZ(cbIdx)-1; z>=0 && (matchFound); z--)
	{
		matchFound = false;

		for(int i=0; i<numIntContsAtZ(z); i++)
		{
			int cidx = idxIntContZ(z,i);
			if( contoursSameSurf(cidx, prevFoundIdx) )
			{
				prevFoundIdx = cidx;                  // set this to "previously found contour"
				matchFound   = true;                  // allows another iteration of the loop
				setDeleteFlag( getC(cidx), 1 );
				numDeleted++;
				//break;
			}
		}
	}

	//## SEARCH UP AND FLAG MATCHING INTERPOLATED CONTOURS FOR DELETION:

	prevFoundIdx = cbIdx;
	matchFound = true;
	for (int z=getCZ(cbIdx)+1; z<interpolation_data.zsize && (matchFound); z++)
	{
		matchFound = false;

		for(int i=0; i<numIntContsAtZ(z); i++)
		{
			int cidx = idxIntContZ(z,i);
			if( contoursSameSurf(cidx, prevFoundIdx) )
			{
				prevFoundIdx = cidx;                  // set this to "previously found contour"
				matchFound   = true;                  // allows another iteration of the loop
				setDeleteFlag( getC(cidx), 1 );
				numDeleted++;
				//break;
			}
		}
	}

	return numDeleted;
}



//-------------------------------
//** SURFACE RESOLVING METHODS:




//------------------------
//-- Returns the two points (one from each contour) which best correspond
//-- by applying the point corrspondance algorithm specified by interpolation_data.ptResolveMethod

void imod_interpolation::
findCoorrespondingPts( Icont *cont1, Icont *cont2,
											 int *idxCont1, int *idxCont2 )
{
	if( interpolation_data.ptResolveMethod == ptmethod::PT_FOUR_PTS )
	{
		findCoorrespondingPts_FourPtMBR    (cont1,cont2,idxCont1,idxCont2);
	}
	else
	{
		findCoorrespondingPts_AllConvexPts (cont1,cont2,idxCont1,idxCont2);
	}
}

//------------------------
//-- Takes two contours and tries to find two points which are MOST
//-- likely to correspond. In this (relatively processing inexpensive)
//-- function it finds the (convex) points in both contours which are
//-- furthest to left, right, top and bottom. It then works out which pair
//-- of points form the closest angles with the center of their the
//-- minimum bounding rectangle around their contour - and returns these
//-- indexes as "coorespoding points".


void imod_interpolation::
findCoorrespondingPts_FourPtMBR( Icont *cont1, Icont *cont2,
																 int *idxCont1, int *idxCont2 )
{
	if( psize( cont1 ) < 3 || psize( cont2 ) < 3 ) {
		idxCont1 = nullptr;
		idxCont2 = nullptr;
		return;
	}

//## GENERATE A BOUNDING BOX AND DETERMINE THE CENTER OF THE BOX FOR BOTH CONTOURS:

	Ipoint cont1ll, cont1ur;
	Ipoint cont2ll, cont2ur;

	imodContourGetBBox( cont1, &cont1ll, &cont1ur );
	imodContourGetBBox( cont2, &cont2ll, &cont2ur );

	Ipoint centerMBRCont1 = line_getPtHalfwayBetween( &cont1ll, &cont1ur );
	Ipoint centerMBRCont2 = line_getPtHalfwayBetween( &cont2ll, &cont2ur );
						// represents the center of the minimum bounding box of cont1 and cont2

//## IDENTIFY WHICH POINTS (INDEXES) CONTAIN THE MIN AND MAX VALUES
//## ALONG X AND Y (CALL THESE "EXTREME POINTS") FOR BOTH CONTOURS:

	Ipoint *cont1pts = imodContourGetPoints( cont1 );
	Ipoint *cont2pts = imodContourGetPoints( cont2 );

	int minXCont1, maxXCont1, minYCont1, maxYCont1;
								// stores the INDEX of point contain the minimum and maximum
								// values for both X and Y in cont1

	for (int i=0; i< psize(cont1);i++) {
		if      ( cont1pts[i].x == cont1ll.x )
			minXCont1 = i;
		else if ( cont1pts[i].x == cont1ur.x )
			maxXCont1 = i;

		if      ( cont1pts[i].y == cont1ll.y )
			minYCont1 = i;
		else if ( cont1pts[i].y == cont1ur.y )
			maxYCont1 = i;
	}

	int minXCont2, maxXCont2, minYCont2, maxYCont2;
									// stores the INDEX of point contain the minimum and maximum
									// values for both X and Y in cont1

	for (int i=0; i<psize( cont2 ); i++) {
		if      ( cont2pts[i].x == cont2ll.x )
			minXCont2 = i;
		else if ( cont2pts[i].x == cont2ur.x )
			maxXCont2 = i;

		if      ( cont2pts[i].y == cont2ll.y )
			minYCont2 = i;
		else if ( cont2pts[i].y == cont2ur.y )
			maxYCont2 = i;
	}


//## CALCULATE THE DIFFERENCE IN ANGLE BETWEEN THE LINES EXENDING FROM
//## THE CENTER OF BOTH CONTOURS TO EXTREME POINTS

	float angleDiffMinX = ABS(line_getAngle2DPos(&cont1pts[minXCont1], &centerMBRCont1)
													- line_getAngle2DPos(&cont2pts[minXCont2], &centerMBRCont2) );
	float angleDiffMaxX = ABS(line_getAngle2DPos(&cont1pts[maxXCont1], &centerMBRCont1)
													- line_getAngle2DPos(&cont2pts[maxXCont2], &centerMBRCont2) );
	float angleDiffMinY = ABS(line_getAngle2DPos(&cont1pts[minYCont1], &centerMBRCont1)
													- line_getAngle2DPos(&cont2pts[minYCont2], &centerMBRCont2) );
	float angleDiffMaxY = ABS(line_getAngle2DPos(&cont1pts[maxYCont1], &centerMBRCont1)
													- line_getAngle2DPos(&cont2pts[maxYCont2], &centerMBRCont2) );


//## DETERMINE WHICH EXTREME POINTS MATCH THE BEST (THOSE WITH THE SMALLEST
//## DIFFERENCE IN ANGLE) AND RETURN THESE "COORESPONDING" POINTS

	float lowestAngleDiff = MIN(angleDiffMaxY,MIN(angleDiffMinY,
																								MIN(angleDiffMaxX,angleDiffMinX)));

	if		(angleDiffMinX == lowestAngleDiff) {
		*idxCont1 = minXCont1;
		*idxCont2 = minXCont2;
	}
	else if	(angleDiffMaxX == lowestAngleDiff) {
		*idxCont1 = maxXCont1;
		*idxCont2 = maxXCont2;
	}
	else if	(angleDiffMinY == lowestAngleDiff) {
		*idxCont1 = minYCont1;
		*idxCont2 = minYCont2;
	}
	else {
		*idxCont1 = maxYCont1;
		*idxCont2 = maxYCont2;
	}
}




//------------------------
//-- Takes two clockwise contours and tries to find two points which are MOST
//-- likely to correspond. For each contour it first finds all the convex
//-- points and determines the angle they make with the center of the MBR.
//-- The list of points indexes and angles are then compared and the indexes
//-- of the two points (one from each contour) with the closest angle
//-- are returned as the as best "coorespoding points".

void imod_interpolation::
findCoorrespondingPts_AllConvexPts( Icont *cont1, Icont *cont2,
																					int *idxCont1, int *idxCont2 )
{
	if( psize( cont1 ) < 3 || psize( cont2 ) < 3 ) {
		idxCont1 = nullptr;
		idxCont2 = nullptr;
		return;
	}

	//## GENERATE A MINIMUM BOUNDING BOX (MBR) AROUND EACH CONTOUR AND DETERMINE IT'S SIZE:

	Ipoint cont1ll, cont1ur;
	Ipoint cont2ll, cont2ur;

	imodContourGetBBox( cont1, &cont1ll, &cont1ur );
	imodContourGetBBox( cont2, &cont2ll, &cont2ur );


	float cont1MBRWidth  = cont1ur.x - cont1ll.x;
	float cont1MBRHeight = cont1ur.y - cont1ll.y;

	float cont2MBRWidth  = cont2ur.x - cont2ll.x;
	float cont2MBRHeight = cont2ur.y - cont2ll.y;

	if( cont1MBRWidth==0 || cont1MBRHeight==0 || cont2MBRWidth==0 || cont2MBRHeight==0 )
	{
		wprint("\aContour with zero area was found");
		findCoorrespondingPts_FourPtMBR( cont1, cont2, idxCont1, idxCont2 );
		return;
	}


	//## CREATE A COPY OF THE CONTOURS AND SET THE Z VALUE OF EACH
	//## POINT TO EQUAL IT'S ORIGINAL INDEX (SO WE DON'T LOSE IT WHEN WE MODIFY):

	Icont *c1 = imodContourDup(cont1);
	Icont *c2 = imodContourDup(cont2);

	for( int p=0; p<psize(c1); p++ )
		getPt(c1,p)->z = p;
	for( int p=0; p<psize(c2); p++ )
		getPt(c2,p)->z = p;


	//## SCALE THE CONTOUR WITH THE SMALLER MINIMUM BOUNDING RECTANGLE (MBR)
	//## SO IT'S MBR BECOMES THE SAME SIZE AS THE LARGER ONE

	float cont1MBRArea   = cont1MBRWidth * cont1MBRHeight;
	float cont2MBRArea   = cont2MBRWidth * cont2MBRHeight;

	if( cont1MBRWidth == cont2MBRWidth &&     // if both MBRs are same dimensions:
			cont1MBRHeight == cont2MBRHeight )    //    do nothing
	{
		;
	}
	if( cont2MBRArea > cont1MBRArea )   // if cont2's MBR is bigger: scale up cont1
	{
		float cont1ScaleX = fDiv( cont2MBRWidth,  cont1MBRWidth );
		float cont1ScaleY = fDiv( cont2MBRHeight, cont1MBRHeight );
		cont_scaleAboutPtXY( c1, &cont1ll, cont1ScaleX, cont1ScaleY );
		imodContourGetBBox( c1, &cont1ll, &cont1ur );

		cont1MBRWidth  = cont1ur.x - cont1ll.x;
		cont1MBRHeight = cont1ur.y - cont1ll.y;
	}
	else                                // else (if cont1's MBR is bigger): scale up cont2
	{
		float cont2ScaleX = fDiv( cont1MBRWidth,  cont2MBRWidth );
		float cont2ScaleY = fDiv( cont1MBRHeight, cont2MBRHeight );
		cont_scaleAboutPtXY( c2, &cont2ll, cont2ScaleX, cont2ScaleY );
		imodContourGetBBox( c2, &cont2ll, &cont2ur );

		cont2MBRWidth  = cont2ur.x - cont2ll.x;
		cont2MBRHeight = cont2ur.y - cont2ll.y;
	}

	Ipoint centerMBRCont1 = line_getPtHalfwayBetween( &cont1ll, &cont1ur );
	Ipoint centerMBRCont2 = line_getPtHalfwayBetween( &cont2ll, &cont2ur );


	//## MAKE BOTH CONTOURS CLOCKWISE AND CONVEX (REMOVE ALL CONCAVE POINTS)

	//imodContourMakeDirection(c1, IMOD_CONTOUR_CLOCKWISE);     //|- make sure conts
	//imodContourMakeDirection(c2, IMOD_CONTOUR_CLOCKWISE);     //|  are clockwise

	cont_makeConvex(c1);
	cont_makeConvex(c2);


	//## FOR EACH (CONVEX) POINT IN BOTH CONTOURS, CALCULATE THE ANGLE IT MAKES
	//## WITH THE CENTER OF THE CONTOUR'S MBR AND POPULATE THIS AS IT'S X VALUE
	//## AND DETERMINE THE POINT WITH THE GREATEST ANGLE:

	int   c1IdxMaxAngle      = -1;
	float maxAngleToCenterC1 = FLOAT_MIN_POS;
	for( int p=0; p<psize(c1); p++ ) {
		float angleTowardsCenter = line_getAngle2DPos( getPt(c1,p), &centerMBRCont1 );
		getPt(c1,p)->x = angleTowardsCenter;
		if( angleTowardsCenter > maxAngleToCenterC1 )
		{
			maxAngleToCenterC1 = angleTowardsCenter;
			c1IdxMaxAngle = p;
		}
	}

	int   c2IdxMaxAngle      = -1;
	float maxAngleToCenterC2 = FLOAT_MIN_POS;
	for( int p=0; p<psize(c2); p++ ) {
		float angleTowardsCenter = line_getAngle2DPos( getPt(c2,p), &centerMBRCont2 );
		getPt(c2,p)->x = angleTowardsCenter;
		if( angleTowardsCenter > maxAngleToCenterC2 )
		{
			maxAngleToCenterC2 = angleTowardsCenter;
			c2IdxMaxAngle = p;
		}
	}


	//## REORDER THE POINTS IN BOTH CONTOURS TO START FROM THE POINT
	//## WITH THE HIGHEST ANGLE AND DUPLICATE THIS ANGLE AT THE END:

	cont_reorderPtsToStartAtIdx( c1, c1IdxMaxAngle );
	cont_reorderPtsToStartAtIdx( c2, c2IdxMaxAngle );

	imodPointAppendXYZ( c1, getFirstPt(c1)->x-360.0f, 0, getFirstPt(c1)->z );
	imodPointAppendXYZ( c2, getFirstPt(c2)->x-360.0f, 0, getFirstPt(c2)->z );


	//## TRAVERSE FROM THE NEW START OF BOTH CONTOURS TO THE END WITH DESCENDING ANGLE
	//## VALUES TO DETERMINE THE TWO POINTS (ONE FROM EACH CONTOUR) WITH THE CLOSEST
	//## ANGLE TO THEIR MBR:

	int cont1IdxMinDiff = -1;
	int cont2IdxMinDiff = -1;

	float minAngleDiff = FLOAT_MAX;

	int p1 = 0;   // the current index in c1
	int p2 = 0;   // the current index in c2


	while( p1 < psize(c1) && p2 < psize(c2) )
	{
		float currAngC1 = getPt(c1,p1)->x;
		float currAngC2 = getPt(c2,p2)->x;
		float angleDiffBetweenPts = ABS( currAngC1 - currAngC2 );

		if( angleDiffBetweenPts < minAngleDiff )
		{
			minAngleDiff    = angleDiffBetweenPts;
			cont1IdxMinDiff = getPt(c1,p1)->z;
			cont2IdxMinDiff = getPt(c2,p2)->z;
		}

		float nextAngInC1 = getPt(c1,p1+1)->x;
		float nextAngInC2 = getPt(c2,p2+1)->x;

		if( nextAngInC1 > nextAngInC2 )   // if next c1 is greater: increpment p1
			p1++;
		else                              // if next c2 is greater: increpment p2
			p2++;
	}

	*idxCont1 = cont1IdxMinDiff;
	*idxCont2 = cont2IdxMinDiff;

	imodContourDelete(c1);
	imodContourDelete(c2);
}





//------------------------
//-- Takes a base contour which branches into two contours and
//-- brakes the base contour into two new "branch" contours which
//-- roughly match to the branched contours, so they can be directly
//-- interpolated/connected together.
//--    _____         _____
//--   /     \       /    _|   - contBranch1 & contBranch2
//--   \____/        \___/
//--     ______   _____
//--    /      \_/     \       - baseCont
//--    \______________/
//--     ______   _____
//--    /      \ /     \       - newBranch1 & newBranch2
//--    \_______|______/

void imod_interpolation::
breakBaseContIntoTwoForBranching( Icont *baseCont,
																			 Icont *contBranch1, Icont *contBranch2,
																			 Icont *newBranch1, Icont *newBranch2 )
{
	Ipoint centerB1;
	cont_getCenterOfMBR(contBranch1, &centerB1);
	Ipoint centerB2;
	cont_getCenterOfMBR(contBranch2, &centerB1);

	float angleSeperatingLine = line_getAngle2DPos(&centerB1,&centerB2) + 90;

	int dummyInt;
	float minDistB1ToCent2;   //|- stores the minimum distance from both branch contour
	float minDistB2ToCent1;		//|  to the centroid of the other one

	Ipoint closestPtBranch1;  //|- stores the closest point in both branch contour
	Ipoint closestPtBranch2;	//|  to the centroid of the other one

	cont_findClosestPtInContToGivenPt( &centerB2, contBranch1,
																		 &minDistB1ToCent2, &closestPtBranch1, &dummyInt );
								// gets the closest point (& distance) in
								// branch 1 contour to the centroid of branch 2
	cont_findClosestPtInContToGivenPt( &centerB1, contBranch2,
																		 &minDistB2ToCent1, &closestPtBranch2, &dummyInt );
								// gets the closest point (& distance) in branch 2
								// contour to the centroid of branch 1

	float distBetweenCentroids = line_distBetweenPts2D( &centerB1, &centerB2 );
	float gapBetweenClosestPts = (minDistB1ToCent2 + minDistB2ToCent1)
																- distBetweenCentroids;
	float fractBetweenCentroidsGapMiddle =
									(minDistB2ToCent1 - (gapBetweenClosestPts/2)) / distBetweenCentroids;
	Ipoint ptMiddleGap = line_findPtFractBetweenPts( &centerB1, &centerB2,
																									fractBetweenCentroidsGapMiddle );

	bool isMiddleGapInsideBase = imodPointInsideCont( baseCont, &ptMiddleGap );
}




//********************************

//## DIFFERENT METHODS FOR TILING TWO KEY CONTOURS FOR INTERPOLATION:

//  The follow functions are different methods to prepare two key
//  contours for interpolation. These functions take two key contours,
//  creates new versions of these such that they both have the same
//  number of points, and these points correspond one-to-one, ready for
//  point interpolation.





//------------------------
//-- Takes two key contours (which it can close if necessary), and
//-- creates new versions of each contour, which have an equal number of
//-- points in corresponding order, so that each point can be easily interpolated.
//--
//-- The algorithm used to determine connecting points depends
//-- on the value of 'interpolation_data.tilingMethod'


void imod_interpolation::
modifyKeyContsForInterp( int clIdx, int cuIdx,
												 Icont *contLNew, Icont *contUNew,
												 bool closeContours,
												 bool findCoorrespondingPtsAndMakeClockwise )
{
	//## APPLY USER SELECTED TILING ALGORITHM:


	switch(interpolation_data.tilingMethod)
	{
		case(tilingmethod::TM_CONSERVE_LENGTH):
		{
			modifyKeyContsForInterp_LineConservation( getC(clIdx), getC(cuIdx),
																								contLNew, contUNew, closeContours,
																								findCoorrespondingPtsAndMakeClockwise );
			break;
		}

		case(tilingmethod::TM_MIN_SA):
		{
			modifyKeyContsForInterp_MinimumAreaCost( clIdx, cuIdx, contLNew, contUNew );
			break;
		}

	}
}


//------------------------
//-- Takes two key contours (which it can close if necessary), and
//-- creates new versions which are ready for one-to-one point interpolation.
//-- It does this adding points to each, such that coorresponding points
//-- on each contour are the same proportional distance along the
//-- length of their contour.
//-- NOTE: This method achieves better effective interpolation
//-- than tiling (mesh fitting) algorithms which do not add extra
//-- points, but instead connect existing ones.

void imod_interpolation::
modifyKeyContsForInterp_LineConservation( Icont *contKeyL, Icont *contKeyU,
																					Icont *contLNew, Icont *contUNew,
																					bool closeContours,
																					bool findCoorrespondingPtsAndMakeClockwise )
{
	imodContourDefault( contLNew );
	imodContourDefault( contUNew );

	if( isEmpty(contKeyL) || isEmpty(contKeyU) || getZ(contKeyL) >= getZ(contKeyU) ) {
		wprint( "\aERROR: modifyKeyContsForInterp_LineConservation()\n" );
		return;
	}

	//## MODIFY KEY CONTOURS SUCH THAT FIRST AND LAST POINTS CORRESPOND

	Icont *contL = imodContourDup(contKeyL);	//|- stores modified version of contour where
	Icont *contU = imodContourDup(contKeyU);  //|  first and last points in each correspond

	imodContourUnique( contL );		// |-- removes any surpurfluous points from contours
	imodContourUnique( contU );		// |

	if( psize(contL) == 1 ) {
		imodContourDelete(contUNew);
		contUNew = imodContourDup( contU );
		for( int i=0; i<psize( contU ); i++ )
			imodPointAppend( contLNew, getPt( contL, 0) );
		imodContourDelete( contL );
		imodContourDelete( contU );
		return;
	}

	if( psize(contU) == 1 ) {
		imodContourDelete(contLNew);
		contLNew = imodContourDup( contL );
		for( int i=0; i<psize( contL ); i++ )
			imodPointAppend ( contUNew, getPt( contU, 0) );
		imodContourDelete( contL );
		imodContourDelete( contU );
		return;
	}

	if(closeContours)        // if we want to close contours:
	{
		imodPointAppend(contL, getPt(contL,0)); // append first point to end
		imodPointAppend(contU, getPt(contU,0));
	}


	if ( findCoorrespondingPtsAndMakeClockwise )
	{
		imodContourMakeDirection(contL, IMOD_CONTOUR_CLOCKWISE);	//|- make sure conts
		imodContourMakeDirection(contU, IMOD_CONTOUR_CLOCKWISE);	//|  are clockwise

		int idxStartContL, idxStartContU;

		findCoorrespondingPts(contL,contU,&idxStartContL,&idxStartContU);
							// finds two points (indexes) on the two contours which appear to match up

		cont_reorderPtsToStartAtIdx( contL, idxStartContL );
		cont_reorderPtsToStartAtIdx( contU, idxStartContU );
	}


	//## CREATE ARRAYS FOR BOTH CONTOURS SHOWING, FOR EACH POINT THE
	//## FRACTION OF ITS DISTANCE ALONG THE TOTAL LENGTH OF THE CONTOUR

	std::vector<float> contLFractsV = cont_getFractPtsAlongLength(contL,false,0);
	std::vector<float> contUFractsV = cont_getFractPtsAlongLength(contU,false,0);


	//## CREATE NEW CONTOURS WITH EQUAL NUMBER OF POINTS WHICH MATCH
	//## UP USING LINE FRACTIONS METHOD:

	cont_addPtsFractsAlongLength( contL, contLNew, contUFractsV, false, true, 0 );
	cont_addPtsFractsAlongLength( contU, contUNew, contLFractsV, false, true, 0 );

	imodContourDelete( contL );
	imodContourDelete( contU );
}




//------------------------
//-- Uses an area minimization metric to creates new versions
//-- of the contours ready for one-to-one point interpolation.
//--
//-- NOTE: It achieves this by first calling David's imodObjectMeshConts function
//-- (see: mkmesh.c > imeshContoursCost) to create a triangular mesh between
//-- the two contours a minimum area cost metric. It then looks through the
//-- lits of vertexs and indexes in the mesh to work out which points in one contour
//-- connect to the points in the other contour and returns two modified version
//-- of the contours which each have an equal number of cooresponding points.

void imod_interpolation::
modifyKeyContsForInterp_MinimumAreaCost( int cbIdx, int ctIdx,
																					Icont *contBNew, Icont *contTNew )
{
	imodContourDefault(contBNew);
	imodContourDefault(contTNew);

	Icont *contB = getC(cbIdx);               // bottom contour
	Icont *contT = getC(ctIdx);               // top contour

	bool inside = 0;

	Ipoint scalePt;
	setPt( &scalePt, 1, 1, 1 );

	Imesh *mesh(nullptr);             // stores a mesh with index connections between
																	// the verticies in contB amd contT
	mesh = imeshContoursCost( obj, contB, contT, &scalePt, inside, cbIdx, ctIdx  );

	int numVert = imodMeshGetMaxVert(mesh);   // number of vertexes
	int numIndx = imodMeshGetMaxIndex(mesh);  // number of indexes


	//## PRINT MESH DATA:
	/*
	wprint( "MESH RESULTS\n");

	wprint("\n  numVert=%d\n", numVert );
	for(int v=0; v<numVert; v++ )
	{
		Ipoint *vert = imodMeshGetVert(mesh, v);
		wprint( "vert %d -> %d,%d,%d \n", v, (int)vert->x, (int)vert->y, (int)vert->z );
	}

	wprint("\n  numIndx=%d\n", numIndx );
	for(int i=0; i<imodMeshGetMaxIndex(mesh); i++ )
	{
		int index = imodMeshGetIndex(mesh, i);
		//wprint( "(%d) %d ,", i, index );
		wprint( "%d ,", index );
	}
	*/


	//## DETERMINE THE FIRST VERTEX IN THE TOP CONTOUR:

	int firstVertT = 0;         // the first vertex in the mesh vertex array belonging
															// to the top contour (all before are the bottom contour)

	for(int v=1; v<numVert; v++ )
		if( imodMeshGetVert(mesh, v-1)->z != imodMeshGetVert(mesh, v)->z )
		{
			firstVertT = v;
			break;
		}

	if( firstVertT != psize(contB)+1 )
		wprint("\aWARNING: modifyKeyContsForInterp_MinimumAreaCost - unexpected # points\n");


	//## MATCH THE VERTEXES TO THE TWO ORIGINAL CONTOURS:
	/*
	vector<IntAndInt> vertMatch;
	vertMatch.resize( numVert );

	for(int v=0; v<firstVertT; v++ )
	{
		Ipoint *vert = imodMeshGetVert(mesh, v);

		for( int p=0; p<psize(contB); p++ )
		{
			if( imodPointIsEqual( vert, getPt(contB,p) ) )
			{
				vertMatch[v].idx1 = 0;
				vertMatch[v].idx2 = p;
			}
		}
	}
	for(int v=firstVertT; v<numVert; v++ )
	{
		Ipoint *vert = imodMeshGetVert(mesh, v);

		for( int p=0; p<psize(contT); p++ )
		{
			if( imodPointIsEqual( vert, getPt(contT,p) ) )
			{
				vertMatch[v].idx1 = 1;
				vertMatch[v].idx2 = p;
			}
		}
	}

	for(int i=0; i<(int)vertMatch.size(); i++ )   //%%%%%%%
		wprint( "vertMatch %d = %d - %d \n", i, vertMatch[i].idx1, vertMatch[i].idx2 );
	*/

	//## CREATE LIST OF CONNECTIONS:

	std::vector<IntAndInt> conn;           // stores a list of connections matching
																		// a point index in contB to a point index in contT

	//int increment = (numVert>1000) ? 6 : 3;
	//int increment = 1;


	for(int i=1; i<imodMeshGetMaxIndex(mesh); i+=1 )
	{
		int vIdx1 = imodMeshGetIndex(mesh, i-1);      // first  vertex index in edge
		int vIdx2 = imodMeshGetIndex(mesh, i);        // second vertex index in edge

		bool eitherIdxInstruction = vIdx1 < 0 || vIdx2 < 0;
		bool bothVIdxContB        = vIdx1  < firstVertT && vIdx2  < firstVertT;
		bool bothVIdxContT        = vIdx1 >= firstVertT && vIdx2 >= firstVertT;

		if ( !eitherIdxInstruction      // if neither vertex index is an instruction
				 && !bothVIdxContB          // and the edge includes one vertex
				 && !bothVIdxContT )        // from each contour
		{
			int idxB = MIN(vIdx1,vIdx2);    // the vertex index in the bottom contour
			int idxT = MAX(vIdx1,vIdx2);    // the vertex index in the top    contour

			conn.push_back( IntAndInt( idxB, idxT ) );   // add to list
		}
	}


	//## SORT LIST OF CONNECTIONS AND ELIMINATE DUPLICATES:

	conn = vector_sort( conn );
	vector_eliminateDuplicates( conn );

	//for(int i=0; i<(int)conn.size(); i++ )   //%%%%%%%
	//  wprint( "connection %d = %d - %d \n", i, conn[i].idx1, conn[i].idx2 );


	//## CREATE MATCHING CONTOURS:

	for(int i=0; i<(int)conn.size(); i+=3 )
	{
		int idxB = conn[i].idx1;
		int idxT = conn[i].idx2;

		if( idxB >= numVert || idxT >= numVert ) {     // should never happen
			wprint("\aERROR: modifyKeyContsForInterp_MinimumAreaCost()\n");
			continue;
		}

		imodPointAppend( contBNew, imodMeshGetVert(mesh,idxB)  );
		imodPointAppend( contTNew, imodMeshGetVert(mesh,idxT)  );
	}

	/*
	for(int p=psize(contBNew); p<psize(contBNew); p--)
	{

	}*/


	//## IF FIRST POINT OF BOTH CONTOURS IS REPEATED AT THE END: REMOVE IT

	if(    imodPointIsEqual( getFirstPt(contBNew), getLastPt(contBNew)  )
			&& imodPointIsEqual( getFirstPt(contTNew), getLastPt(contTNew)  ) )
	{
		imodPointDelete( contBNew, psize(contBNew)-1 );
		imodPointDelete( contTNew, psize(contTNew)-1 );
	}

	conn.clear();

	imodMeshDelete(mesh);
}







//------------------------
//-- Takes an ordered sequence of four contours and generates "crude" smoothly
//-- interpolated contours between c1 and c2. The interpolation is crude
//-- in that contours themselves are interpolated linearly, but then
//-- their centroid and size of each interpolated contour is adjusted
//-- using catmull rom-spline algorithm over the four key contour.


void imod_interpolation::
interp_SmoothCrude_BetweenConts( int c0, int c1, int c2, int c3 )
{
	Icont *cont0 = getC( c0 );
	Icont *cont1 = getC( c1 );
	Icont *cont2 = getC( c2 );
	Icont *cont3 = getC( c3 );

	int c1Z = getZInt( cont1 );
	int c2Z = getZInt( cont2 );
	int diffInSlices = c2Z - c1Z;				// the difference in Z between the two key contours

	if( diffInSlices <=1 )
		return;


	//## FOR EACH KEY CONTOUR: CALCULATE AND STORE IT'S CENTROID AND RADIUS:

	Ipoint c0Centroid, c1Centroid, c2Centroid, c3Centroid;

	cont_getCentroid( cont0, &c0Centroid );
	cont_getCentroid( cont1, &c1Centroid );
	cont_getCentroid( cont2, &c2Centroid );
	cont_getCentroid( cont3, &c3Centroid );

	float c0Radius = cont_getRadius( cont0 );
	float c1Radius = cont_getRadius( cont1 );
	float c2Radius = cont_getRadius( cont2 );
	float c3Radius = cont_getRadius( cont3 );


	//## GENERAE A SET OF LINEARLY INTERPOLATED CONTOURS BETWEEN CONT1 AND CONT2:

	std::vector<IcontPtr> newConts = getLinearInterpConts( c1, c2 );

	std::vector<float> fractsZ = calcCardinalSplineFractsEachSlice( c0Centroid, c1Centroid,
																																	c2Centroid, c3Centroid );
					// stores a vector of fractions representing how far between
					// cont1 and cont2 a cardinal spline between the centroids crosses
					// EACH z slice between them (in ascending order)


	//## FOR EACH INTERPOALTED CONTOUR MOVE AND RESIZE IT USING THE CARDINAL SPLINE
	//## VALUES FOR CENTROID AND RADIUS:

	for (int c=0; c<int(newConts.size()); c++)			 // for each interpolated contour:
	{
		Icont *cont = newConts[c].cont;
		int interpCurrZ = getZInt(cont);

		float fractUpFromFirst = fDiv( float(interpCurrZ - c1Z), float(diffInSlices) );
						// the fraction distance the current Z value
						// is between lower and upper key contours

		//## TRANSLATE CONTOUR USING CARDINAL SPLINE THROUGH CENTROID:

		Ipoint currCentroid;
		//cont_getCentroid( cont, &currCentroid );
		currCentroid = line_findPtFractBetweenPts(&c1Centroid,&c2Centroid,fractUpFromFirst);
												// the current centroid of the interpolated contour

		Ipoint newCentroid = getPtCardinalSpline( fractsZ[c], c0Centroid, c1Centroid,
																							c2Centroid, c3Centroid, TENSILE_FRACT );
											// the centroid we want (according to the spline algorithm)

		cont_translate( cont, newCentroid.x-currCentroid.x, newCentroid.y-currCentroid.y );
											// translates (moves) the interpolated contour
											// so the centoid is in the correct position

		//## SCALE CONTOUR USING CARDINAL SPLINE THROUGH RADIUS:

		float currRadius = cont_getRadius( cont );

		float newRadius = getValCardinalSpline( fractsZ[c], c0Radius, c1Radius,
																						c2Radius, c3Radius, TENSILE_FRACT );
											// the radius we want (according to the spline algorithm)

		float scaleFactor = fDiv( newRadius, currRadius );
											// factor by which we will resize the interpolated contour

		cont_scaleAboutPt( cont, &newCentroid, scaleFactor, true );
											// resizes the contour about it's centroid to
											// the desired size (represented by radius)


		//## ADD THE INTERPOLATED CONTOUR TO THE OBJECT:

		setZValue( cont, interpCurrZ );
		addInterpolatedContToObj( interpolation_data, obj, cont );
	}

	deleteContours( newConts );
}






//------------------------
//-- Takes two key contours ("contKeyL", "contKeyU") and
//-- generates a series of smooth interpolated contours between them.
//-- It does this by trying to find contour c0 below contKeyL
//-- and cont c3 above contKeyU, then making versions of these
//-- such that they all have the same # of points, then interpolating
//-- between each set of cooresponing points using cardinal splines.

void imod_interpolation::
interp_SmoothPointwise_BetweenTwoConts( int clIdx, int cuIdx )
{
	Icont* contL = getC(clIdx);
	Icont* contU = getC(cuIdx);

	if( psize(contL) == 0 || psize(contU) == 0
			|| getZ(contL) >= getZ(contU) ) {
		std::cerr << "ERROR: interp_SmoothPointwise()" << std::endl;
		return;
	}
	int zDistApart = getZ(contU) - getZ(contL);
	if( zDistApart <= 1 )				// if conts only one slice apart, or lower contour
		return;                   // is above upper contour: return

//## PREPARE UPPER AND LOWER KEY CONTOUR FOR LINEAR INTERPOLATION BY
//## ADDING EXTRA POINTS TO EACH:

	Icont *c0 = imodContourNew();     //|- will store final contour, all with an
	Icont *c1 = imodContourNew();     //|  equal number of points, ready for
	Icont *c2 = imodContourNew();     //|  smooth interpolation.
	Icont *c3 = imodContourNew();     //|

	modifyKeyContsForInterp( clIdx, cuIdx, c1, c2, true, true );
							// prepares the lower and upper key contour (now called c1 and c2)
							// for linear interpolation by adding extra points to each
							// such that they an equal number of (cooresponding) points

//## IF THERE IS A CONTOUR BELOW THE (NOW MODIFIED) LOWER CONTOUR:
//## CREATE A VERSION WITH COORESPONDING POINTS READY FOR SMOOTH INTERPOLATION

	int idxC0 = findIdxNearestKeyCont( clIdx, interpolation_data.zBridge, false );
	if(idxC0 == -1)     // if no contour was found below c1: c0 is just c1
	{
		imodContourDelete(c0);
		c0 = imodContourDup(c1);
	}
	else                // else: generate a version of this into c0 which has
	{                   // points cooresponding one-to-one with c2
		int idxStartC0, idxStartC1;

		findCoorrespondingPts( getC(idxC0), c1,&idxStartC0,&idxStartC1);
		std::vector<float> c1FractsV = cont_getFractPtsAlongLength(c1,true,idxStartC1);
		cont_addPtsFractsAlongLength( getC(idxC0), c0, c1FractsV, true,false,idxStartC0);
		cont_reorderPtsToStartAtIdx( c0, psize(c1)-idxStartC1 );
	}


//## IF THERE IS A CONTOUR ABOVE THE (NOW MODIFIED) UPPER CONTOUR:
//## CREATE A VERSION WITH COORESPONDING POINTS READY FOR SMOOTH INTERPOLATION

	int idxC3 = findIdxNearestKeyCont( cuIdx, interpolation_data.zBridge, true );
	if(idxC3 == -1)		// if no contour was found above c2: c3 is just c2
	{
		imodContourDelete(c3);
		c3 = imodContourDup(c2);
	}
	else            // else (if contour was found above c2): generate a version of
	{               // this into c3 which has points cooresponding one-to-one with c2
		int idxStartC2, idxStartC3;

		findCoorrespondingPts_FourPtMBR   ( c2, getC(idxC3), &idxStartC2, &idxStartC3 );
		std::vector<float> c2FractsV = cont_getFractPtsAlongLength(c2,true,idxStartC2);
		cont_addPtsFractsAlongLength(getC(idxC3), c3, c2FractsV, true,false,idxStartC3);
		cont_reorderPtsToStartAtIdx( c3, psize(c2)-idxStartC2 );
	}


	if( (psize(c0) != psize(c1)) || (psize(c1) != psize(c2)) || (psize(c2) != psize(c3)) )
	{
		wprint("\aERROR: interp_SmoothPointwise() - conts not all same # points\n");
		//cout << "psize(c0)=" << psize(c0) << endl;     flush(std::cout); //%%%%%%%%%
		//cout << "psize(c1)=" << psize(c1) << endl;     flush(std::cout); //%%%%%%%%%
		//cout << "psize(c2)=" << psize(c2) << endl;     flush(std::cout); //%%%%%%%%%
		//cout << "psize(c3)=" << psize(c3) << endl;     flush(std::cout); //%%%%%%%%%
		return;
	}


//## USING THE FOUR COORESPONDING COUNTOURS (c0,c1,c2,c3):
//## INTERPOLATE EACH COORESPONDING POINT BETWEEN c1 AND c2
//## TO FORM "SMOOTH" INTERPOLATED CONTOURS

	std::vector<IcontPtr> newConts;				// will store the new (interpolated) contours
																		// between the upper and lower contours (c1 and c2)

	for(int i=0; i<(zDistApart-1); i++  )	// populate newConts with empty contours
		newConts.push_back( IcontPtr() );

	for (int j=0; j<psize(c1); j++)		// for each (cooresponding) point:
	{
		std::vector<float> fractsZ =
			calcCardinalSplineFractsEachSlice( *getPt(c0,j), *getPt(c1,j),
																											 *getPt(c2,j), *getPt(c3,j)  );
						// fractsZ stores a vector of fractions representing how far
						// between the lower and upper key contours  a cardinal spline between
						// the centroids crosses EACH z slice between them (in ascending order)

		for(int i=0; i<int(newConts.size()); i++)	// for each z between upper & lower slice:
		{
			Ipoint newPt = getPtCardinalSpline( fractsZ[i],  *getPt(c0,j), *getPt(c1,j),
																			*getPt(c2,j), *getPt(c3,j), TENSILE_FRACT  );
														// finds the point at which the point intersects the
														// desired z value (according to the catmull-rom algorithm)

			newPt.z = getZ(c1)+1 + i;
														// changes the z value so it's the
														// correct z slice (it might be slightly off)

			imodPointAppend( newConts[i].cont, &newPt );
														// add this point to the correct contour.
		}
	}

	for(int i=0; i<int(newConts.size()); i++  )        // for each newly generated contour:
		addInterpolatedContToObj( interpolation_data, obj, newConts[i].cont );	// add it to the object

	deleteContours(newConts);
	imodContourDelete( c0 );
	imodContourDelete( c1 );
	imodContourDelete( c2 );
	imodContourDelete( c3 );
}

//------------------------
//-- Takes two key contours and generates a new set of linearly interpolated
//-- contours between them.

std::vector<IcontPtr> imod_interpolation::
getLinearInterpConts( int c1Idx, int c2Idx )
{
	std::vector<IcontPtr> newInterpolatedConts;

	bool closed = getClosed(c1Idx) || getClosed(c2Idx);
										// typically both (or at least one) of the contours will be
										// closed (used in modifyKeyContsForInterp)

	bool c1IsLower = getCZ(c1Idx) < getCZ(c2Idx);
	int clIdx     = (c1IsLower) ? c1Idx : c2Idx;
	int cuIdx     = (c1IsLower) ? c2Idx : c1Idx;
	Icont *contL  = getC(clIdx);
	Icont *contU  = getC(cuIdx);


	int zDistApart = getZ(contU) - getZ(contL);
	if( zDistApart <= 1 )   // if conts only one slice apart, or lower is above upper:
	{
		//wprint( "ERROR: getLinearInterpConts()\n" );      // should avoid this
		return newInterpolatedConts;                      // return empty set
	}

//## CREATE TWO NEW KEY CONTOURS WITH (AN EQUAL NUMBER OF) COORESPONDING
//## POINTS FOR INTERPOLATION:

	Icont *contLNew = imodContourNew();		// |-- version of conts with an equal number of
	Icont *contUNew = imodContourNew();		// |   points which coorspond one-to-one.

	modifyKeyContsForInterp( clIdx, cuIdx, contLNew, contUNew, closed, closed );


//## USE NEW KEY CONTOURS TO PERFORM LINEAR POINT INTERPOLATION AND GENERATE
//## NEW INTERPOLATED CONTOURS BETWEEN THE KEY CONTOURS:


	for(int z=getZ( contLNew )+1; z<=getZ( contUNew )-1; z++  )
	{
		float fractBetweenKeyContour = float(z - getZ( contLNew )) / float(zDistApart);

		Icont *contNew = imodContourNew();
		setInterpolated( contNew, 1 );
		for(int i=0; i<psize( contLNew ) && i<psize( contUNew ); i++) {
			Ipoint tmpPt = line_findPtFractBetweenPts( getPt(contLNew,i), getPt(contUNew,i),
									fractBetweenKeyContour );
			imodPointAppend( contNew, &tmpPt);
		}

		newInterpolatedConts.push_back( IcontPtr(contNew) );
		imodContourDelete( contNew );
	}

	imodContourDelete( contLNew );
	imodContourDelete( contUNew );

	if( !closed && isOpenFlag( getC(c1Idx) ) && isOpenFlag( getC(c2Idx) )  )
	{
		//wprint("TWO OPEN CONTS\n");
		for (int i=0; i<newInterpolatedConts.size(); i++)
			setOpenFlag( newInterpolatedConts[i].cont, 1 );
	}


	return newInterpolatedConts;
}


//------------------------
//-- Takes two key contours and generates a new set of linearly
//-- interpolated contours between them.

int imod_interpolation::
interp_Linear_BetweenTwoConts( int clIdx, int cuIdx )
{
	std::vector<IcontPtr> newConts = getLinearInterpConts(clIdx, cuIdx);
	int numNewConts = newConts.size();

	for(int i=0; i<int(newConts.size()); i++)
		addInterpolatedContToObj(interpolation_data, obj, newConts[i].cont );

	deleteContours( newConts );
	return ( numNewConts );
}




//------------------------
//-- Takes a number of key contours and generates a set of linearly
//-- interpolated contours using branching
//-- Returns the number of contours that result

int imod_interpolation::
mergeAllTouchingConts( std::vector<IcontPtr> conts )
{
	for(int i=1; i<int(conts.size()); i++)
		if( cont_doContsTouch( conts[i-1].cont, conts[i].cont ) )
		{
			cont_getOuterUnionPolygon( conts[i-1].cont, conts[i-1].cont, conts[i].cont, 0 );
			deleteAllPts( conts[i].cont );
			eraseContour( conts, i );
			i=0;
		}

	return ( conts.size() );
}

//------------------------
//-- Adds interpolated contours and merges any contours
//-- which are on the same Z slice an touch

void imod_interpolation::
addAllInterpolatedConstAndMerge( std::vector< std::vector<IcontPtr> > newConts )
{
	vector2D_transpose( newConts );

	for (int z=0; z<(int)newConts.size(); z++) {
		int numConts = mergeAllTouchingConts( newConts[z] );
		newConts[z].resize( numConts );   // this line shouldn't be needed, but it is.  :-/
	}

	for (int z=0; z<(int)newConts.size(); z++) {
		for (int b=0; b<(int)newConts[z].size(); b++)
			addInterpolatedContToObj( interpolation_data, obj, newConts[z][b].cont );

		deleteContours( newConts[z] );
	}
}

//------------------------
//-- Takes a number of key contours and generates a set
//-- of linearly interpolated contours using branching

void imod_interpolation::
interp_Linear_BetweenContsMerge( int cbIdx, std::vector<int> branchContIdx )
{
	if( (int)branchContIdx.size() == 1 ) {
		interp_Linear_BetweenTwoConts( cbIdx, branchContIdx[0] );
		return;
	}

	std::vector< std::vector<IcontPtr> > newConts;

	for (int b=0; b<(int)branchContIdx.size(); b++)
		newConts.push_back( getLinearInterpConts(cbIdx, branchContIdx[b]) );

	addAllInterpolatedConstAndMerge( newConts );
}


//------------------------
//-- Takes a single contour, then finds all finds ALL key contours
//-- which (appear to) belong to the same surface and generates
//-- a CRUDE form of smooth interpolation between each of these.
//-- The interpolation is crude in that contours themselves are
//-- interpolated linearly, but then their centroid and size is
//-- adjusted linearly using catmull rom-spline algorithm.
//-- This can be effective for simple contour shapes, but results
//-- are not ideal for more complex-shaped contours.

int imod_interpolation::
interp_Linear( int baseContIdx, int maxDist )
{
	int numAdded = 0;

	switch( interpolation_data.branchingMethod )
	{
		case( branchingmethod::BR_BRANCHING_OFF ):
			{
				if( isKeyContAbove() )
					numAdded += interp_Linear_BetweenTwoConts(sContIdx,closestAboveIdx);
				if( isKeyContBelow() )
					numAdded += interp_Linear_BetweenTwoConts(closestBelowIdx,sContIdx);
			}
			break;


		case( branchingmethod::BR_MERGE_CONTS ):
			{
				if( isKeyContAbove() )
					interp_Linear_BetweenContsMerge( sContIdx, aboveIdxs );
				if( isKeyContBelow() )
					interp_Linear_BetweenContsMerge( sContIdx, belowIdxs );
			}
			break;

		/*
		case( branchingmethod::BR_BRIDGE_GAPS ):
		{
			if( isKeyContAbove() )
			{
				Icont joinedCont = cont_joinContsAtClosestApproach(getAboveConts(obj),true);
				interp_Linear_BetweenTwoConts( obj, getCon(obj), joinedCont );
			}
			if( isKeyContBelow() )
			{
				Icont joinedCont = cont_joinContsAtClosestApproach(getBelowConts(obj),true);
				interp_Linear_BetweenTwoConts( obj, getCon(obj), joinedCont );
			}

		}
		break;

		case( BR_PT_BELOW ):
		{
			wprint("WARNING: BR_PT_BELOW does not work for LINEAR interpolation");
		}
		break;
		*/

		default:
			wprint("Change branching type\n");
			break;
	}

	return numAdded;
}








//------------------------
//-- Takes just one key contour - the "middle contour" in an
//-- approximately spherical object - and generates interpolated contours
//-- above and below it in a sphere shape.


void imod_interpolation::
interp_Spherical_OnSingleCont( int cMiddleIdx, const float modScaleZ )
{
	//Imod *imod = ivwGetModel(interpolation_data.view);
	//float modScaleZ = imodGetZScale(imod);

	Icont *cont = imodContourDup( getC( cMiddleIdx ) );

//## CALCULATE CENTROID, AREA, AND SPAN FOR SPHERE:

	int contZ = getZ( cont );



	float area = imodContourArea(cont);			// will store the area of the middle contour
	Ipoint centroidPt;
	cont_getCenterOfMBR( cont, &centroidPt );
	float radius = sqrt( ABS(area) / PI);			// represents the "average radius" of the
																						// "sphere of best fit" over the contour
	int    zSlicesRadiusSpans = (int)floor(radius/modScaleZ);
																					// number of whole z slices the radius spans


//## DELETE INTERPOLATED CONTOURS IN SPHERE'S WAY:

	deleteAllSameSurfInterpContsInZRange( contZ-(zSlicesRadiusSpans+2),
																				contZ+(zSlicesRadiusSpans+2), cMiddleIdx );


//## GENERATE INTERPOLATED CONTOURS ABOVE AND BELOW MIDDLE CONTOUR:

	for(int i=1; i<=zSlicesRadiusSpans; i++)
	{
		float radiusAtThisZ = sqrt( SQ(radius) - SQ(i*modScaleZ) );
														// using pythagorean theorem:     adj = sqrt(hpy^2 - opp^2)
		float scaleFactor = radiusAtThisZ/radius;
																	// the size of the radius/area formed if the sphere
																	// is cut at this z compared to the sphere's radius
		Icont *contNewLower = imodContourDup( cont );
		cont_scaleAboutPt( contNewLower, &centroidPt, scaleFactor, true );
																	// duplicate middle contour, but resize it about
																	// its centroid so it has the correct radius to
																	// fit the sphere's profile
		Icont *contNewUpper = imodContourDup( contNewLower );

		setZValue(contNewUpper,contZ+i);           //|- add this same scaled down
		addInterpolatedContToObj(interpolation_data, obj, contNewUpper);		//|  contour appropriate distance
		setZValue(contNewLower,contZ-i);           //|  above and below the middle
		addInterpolatedContToObj(interpolation_data, obj, contNewLower);		//|  contour

		imodContourDelete( contNewLower );
		imodContourDelete( contNewUpper );
	}

	imodContourDelete( cont );
}






//------------------------
//-- Takes a single contour, and executes spherical interpolation
//-- by finding ALL key contours which (seem to) belong to the same
//-- surface, determining a "sphere of best fit", and generating
//-- interpolated contours between all pairs of contours, as well as
//-- additional contours above and below the top-most and bottom most
//-- contours respectively (in order to "cap" the sphere).

void imod_interpolation::
interp_Spherical( int baseContIdx, int maxDist, const float modScaleZ )
{
	//Imod *imod = ivwGetModel(interpolation_data.view);
	//float modScaleZ = imodGetZScale(imod);

	if ( interpolation_data.surfResolveMethod != surfacemethod::SR_CENTER_OVERLAP ) {
		wprint( "\r\nNOTE: Am changing surface resolution method to "
						"'centroid overlap' for spherical interpolation.\n" );
		interpolation_data.surfResolveMethod = surfacemethod::SR_CENTER_OVERLAP;
	}

	std::vector<int> idxInSurf =
		findIdxsAllKeyContsNoBranching(baseContIdx,maxDist);

//## IF THERE ARE NO OTHER KEY CONTOURS BELONGING TO THIS SURFACE:
//## ASSUME IT'S THE MIDDLE CONTOUR AND DRAW SPHERE AROUND IT

	if( idxInSurf.size() == 1 ) {
		interp_Spherical_OnSingleCont( baseContIdx, modScaleZ );
		return;
	}


//## (IF THERE WERE OTHER KEY CONTOURS BELONGING TO THIS SURFACE):
//## DELETE ALL INTERPOLATED CONTOURS ATTCHED TO THESE KEY CONTOUR

	for (int i=0; i<int(idxInSurf.size()); i++)
		deleteAllSameSurfInterpContsEitherSideOfCont( idxInSurf[i] );

//## FOR EACH KEY CONTOUR: CALCULATE AND STORE IT'S CENTROID AND RADIUS

	std::vector<Ipoint> contCenter;	// store centroid for each key contour in surface
	std::vector<float>  contRadius;		// store radius   for each key contour in surface

	for (int i=0; i<(int)idxInSurf.size(); i++)		// for each contour idx in same surface:
	{
		Ipoint centroid;                      // center of contour
		float area = imodContourArea( getC( idxInSurf[i] ) );
		cont_getCenterOfMBR( getC(idxInSurf[i] ), &centroid );
		float radius = sqrt(ABS(area)/PI);		// calculate radius for circle with given area

		contCenter.push_back( centroid );     //|- store these values (for this key contour)
		contRadius.push_back( radius );				//|  in an array
	}


//## FOR EACH PAIR OF ADJACENT CONTOURS: CALCULATE AND STORE A CENTER
//## AND RADIUS DESCRIBING A SPHERE ENCOMPASSING THEM

	std::vector<float> sphereZCentroid;		// stores the centroid for a sphere encompassing
																		// the key contour in idxInSurf and the one above it
	std::vector<float> sphereRadius;       // stores the radius for the above

	for (int i=0; i<(int(idxInSurf.size())-1); i++)
	{
		float diffRadius = contRadius[i+1] - contRadius[i];
		float diffZ      = contCenter[i+1].z - contCenter[i].z;

		float radiusMid	= (contRadius[i+1] + contRadius[i]) / 2.0;
		float zMid		= (contCenter[i+1].z + contCenter[i].z) / 2.0;

		float zCenter = (((radiusMid/(diffZ*modScaleZ)) * diffRadius)/modScaleZ) + zMid;
		float zDist   = (contCenter[i].z-zCenter)*modScaleZ;
		float sphereRad = sqrt( SQ(contRadius[i]) + SQ(zDist) );

		sphereZCentroid.push_back ( zCenter   );
		sphereRadius.push_back    ( sphereRad );
	}


//## FOR EACH PAIR OF ADJACENT CONTOURS: GENERATE LINEARLY INTERPOLATED
//## CONTOURS BETWEEN THEM, AND RESIZE THEM USING SPHERE RADIUS

	std::vector<IcontPtr> newConts;

	for (int i=0; i<(int(idxInSurf.size())-1); i++)
	{
		int cIdxStart = idxInSurf[i];
		int cIdxEnd   = idxInSurf[i+1];

		Icont *contStart = getC( idxInSurf[i] );
		Icont *contEnd   = getC( idxInSurf[i+1] );

		newConts = getLinearInterpConts( cIdxStart, cIdxEnd );

		float startZ = getZ(contStart);
		float endZ   = getZ(contEnd);
		float diffInSlices = endZ - startZ;

		//cout << "cIdxStart=" << cIdxStart << " cIdxEnd=" << cIdxEnd << endl;
		//cout << "newConts.size()=" << newConts.size() << " diffInSlices=" << diffInSlices << endl;
		//cout << "startZ=" << startZ << " endZ=" << endZ << endl;


		for (int j=0; j<(int)newConts.size(); j++)
		{
			float interpCurrZ = getZ(newConts[j].cont);
			float fractUpFromFirst = fDiv( (interpCurrZ - startZ), diffInSlices );
			float currRadius = (contRadius[i+1]-contRadius[i])*fractUpFromFirst+contRadius[i];
			float zDist = (sphereZCentroid[i] - interpCurrZ)*modScaleZ;
			float newRadius = sqrt( SQ(sphereRadius[i]) - SQ(zDist) );
														// using pythagorean theorem:     adj = sqrt(hpy^2 - opp^2)
			float scaleFactor = fDiv( newRadius, currRadius );
			Ipoint centroid = line_findPtFractBetweenPts( &contCenter[i], &contCenter[i+1],
																										fractUpFromFirst );
			cont_scaleAboutPt( newConts[j].cont, &centroid, scaleFactor, true );

			addInterpolatedContToObj( interpolation_data, obj, newConts[j].cont );
		}

		deleteContours(newConts);
	}


//## USE SAME METOHD TO CALCULATE A SPHERE FOR THE HIGHEST AND
//## LOWEST CONTOUR (CALL THIS HIGH-LOW-SPHERE):

	int idxLastCont = int( idxInSurf.size() )-1;

	Icont *lowestCont  = getC( idxInSurf[0] );
	Icont *highestCont = getC( idxInSurf[idxLastCont] );

	float diffRadius = contRadius[idxLastCont] - contRadius[0];
	float diffZ      = contCenter[idxLastCont].z - contCenter[0].z;

	float radiusMid	= (contRadius[idxLastCont]-contRadius[0] )/2.0f + contRadius[0];
	float zMid		  = (contCenter[idxLastCont].z-contCenter[0].z)/2.0f + contCenter[0].z;

	float HLsphereZCenter = (((radiusMid/(diffZ*modScaleZ))*diffRadius)/modScaleZ) + zMid;
														// stores the centroid of the "high-low sphere"
														// (a sphere encompassing the top and bottom-most conts)

	float zDist = (contCenter[0].z-HLsphereZCenter)*modScaleZ;
	float HLsphereRadius = sqrt( SQ(contRadius[0]) + SQ(zDist) );		// radius

	float scaleZTEST = ABS(pow((-1.0 * radiusMid * diffRadius /
															((zMid - HLsphereZCenter) * diffZ)), (0.5)));


//## USE HIGH-LOW-SPHERE TO INTERPOLATE UPWARDS FROM HIGHEST AND
//## DOWNWARDS FROM LOWEST CONTOURS:

	int lowestZInt  = (int)ceil (HLsphereZCenter - (HLsphereRadius/modScaleZ));
	int highestZInt = (int)floor(HLsphereZCenter + (HLsphereRadius/modScaleZ));
															// represents the highest and lowest slic where an
															// interpolated contour will be generated


	for(int z=lowestZInt; z<getZ(lowestCont); z++ )				// for each slice below the
	{                                                     // lowest key contour:
		float opp = ((HLsphereZCenter - z)*modScaleZ);
		float newRadius = sqrt( (HLsphereRadius*HLsphereRadius) - (opp*opp) );
		float scaleFactor = fDiv( newRadius, contRadius[0] );
						// determines how small interpolated contour will be compared to key contour

		newConts.push_back( IcontPtr( getC( idxInSurf[0] ) ) );
		cont_scaleAboutPt(newConts.back().cont, &contCenter[0], scaleFactor, true );
											// copy lowest contour, but scale it down to appropriate size ...
		setZValue( newConts.back().cont, z );      // then copy it to desired z slice
	}

	for(int z=highestZInt; z>getZ(highestCont); z-- )			// for each slice above the
	{                                                     // lowest key contour:
		float opp = ((HLsphereZCenter - z)*modScaleZ);
		float newRadius = sqrt( (HLsphereRadius*HLsphereRadius) - (opp*opp) );
		float scaleFactor = fDiv( newRadius, contRadius[idxLastCont] );
						// determines how small interpolated contour will be compared to key contour

		newConts.push_back( IcontPtr( getC( idxInSurf[idxLastCont] ) ) );
		cont_scaleAboutPt(newConts.back().cont,&contCenter[idxLastCont],scaleFactor,true);
											// copy lowest contour, but scale it down to appropriate size ...
		setZValue( newConts.back().cont, z );     // then copy it to desired z slice
	}

	for(int i=0; i<int(newConts.size()); i++)
		addInterpolatedContToObj( interpolation_data, obj, newConts[i].cont );

	deleteContours(newConts);
}







//------------------------
//-- Takes a single contour, then tries to find three key contours above and below
//-- which (appear to) belong to the same surface and generates a CRUDE form of
//-- smooth interpolation between each of these. The interpolation is crude
//-- in that contours themselves are interpolated linearly, but then
//-- their centroid and size is adjusted linearly using catmull rom-spline
//-- algorithm. This can be effective for simple contour shapes, but results
//-- are not ideal for more complex-shaped contours.

void imod_interpolation::
interp_SmoothCrude( int baseContIdx, int maxDist )
{
	if( interpolation_data.branchingMethod == branchingmethod::BR_MERGE_CONTS &&
			( (int)aboveIdxs.size() > 1 || (int)belowIdxs.size() > 1 ) )
	{
		if( isKeyContAbove() )
			interp_Linear_BetweenContsMerge( sContIdx, aboveIdxs );
		if( isKeyContBelow() )
			interp_Linear_BetweenContsMerge( sContIdx, belowIdxs );
		wprint("LINEAR branching instead\n");
		return;
	}

	if( isKeyContAbove() )
	{
		int closestAbove2Idx = findIdxNearestKeyCont(closestAboveIdx,
																								 interpolation_data.zBridge,true);

		int c0 = ( isKeyContBelow() ) ? closestBelowIdx : baseContIdx;
		int c1 = baseContIdx;
		int c2 = closestAboveIdx;
		int c3 = ( closestAbove2Idx!=-1 ) ? closestAbove2Idx : closestAboveIdx;

		interp_SmoothCrude_BetweenConts( c0, c1, c2, c3 );
	}

	if( isKeyContBelow() )
	{
		int closestBelow2Idx = findIdxNearestKeyCont(closestBelowIdx,
																								 interpolation_data.zBridge,false);

		int c0 = ( closestBelow2Idx!=-1 ) ? closestBelow2Idx : closestBelowIdx;
		int c1 = closestBelowIdx;
		int c2 = baseContIdx;
		int c3 = ( isKeyContAbove() ) ? closestAboveIdx : baseContIdx;

		interp_SmoothCrude_BetweenConts( c0, c1, c2, c3 );
	}

}



//------------------------
//-- Takes a single contour, then finds all finds ALL key contours which
//-- (appear to) belong to the same surface and generates a CRUDE form of
//-- smooth interpolation between each of these. The interpolation is crude
//-- in that contours themselves are interpolated linearly, but then
//-- their centroid and size is adjusted linearly using catmull rom-spline
//-- algorithm. This can be effective for simple contour shapes, but results
//-- are not ideal for more complex-shaped contours.

void imod_interpolation::
interp_SmoothCrudeAll( int baseContIdx, int maxDist )
{
	if( interpolation_data.branchingMethod == branchingmethod::BR_MERGE_CONTS &&
			( (int)aboveIdxs.size() > 1 || (int)belowIdxs.size() > 1 ) )
	{
		if( isKeyContAbove() )
			interp_Linear_BetweenContsMerge( sContIdx, aboveIdxs );
		if( isKeyContBelow() )
			interp_Linear_BetweenContsMerge( sContIdx, belowIdxs );
		wprint("LINEAR branching instead\n");
		return;
	}

	std::vector<int> idxInSurf =
			findIdxsAllKeyContsNoBranching( baseContIdx, maxDist );

//## IF THERE ARE NO OTHER KEY CONTOURS BELONGING TO THIS SURFACE: DO NOTHING

	if( idxInSurf.size() == 1 ) {
		return;
	}

//## DELETE ALL INTERPOLATED CONTOURS EITHER SIDE OF THIS KEY CONTOUR

	for (int i=0; i<int(idxInSurf.size()); i++)
		deleteAllSameSurfInterpContsEitherSideOfCont( idxInSurf[i] );


//## FOR EACH PAIR OF ADJACENT CONTOURS: GENERATE CRUDE SMOOTHLY INTERPOLATED
//## CONTOURS BETWEEN THEM

	int numKeyConts = (int)idxInSurf.size();

	for (int i=0; i<=numKeyConts-2; i++)
	{
		int c0	= (i<=0) ? idxInSurf[0] : idxInSurf[i-1];
		int c1	= idxInSurf[i];
		int c2	= idxInSurf[i+1];
		int c3	= (i>=numKeyConts-2) ? idxInSurf[i+1] : idxInSurf[i+2];

		interp_SmoothCrude_BetweenConts( c0, c1, c2, c3 );
	}
}



//------------------------
//-- Takes a single contour, checks if there are other contours
//-- in the same surface, and if so, generates smooth interpolated
//-- contours between them using cardinal splines between
//-- cooresponding points over the contours.

void imod_interpolation::
interp_SmoothPointwise( int baseContIdx, int maxDist )
{
	if( isKeyContAbove() )	{
		interp_SmoothPointwise_BetweenTwoConts(sContIdx,closestAboveIdx);

		int closestAbove2Idx = findIdxNearestKeyCont(closestAboveIdx, interpolation_data.zBridge, true);
		if( closestAbove2Idx!=-1 )
			interp_SmoothPointwise_BetweenTwoConts(closestAboveIdx, closestAbove2Idx);
	}
	if(  isKeyContBelow()  )	{
		interp_SmoothPointwise_BetweenTwoConts(closestBelowIdx,sContIdx);

		int closestBelow2Idx = findIdxNearestKeyCont(closestBelowIdx, interpolation_data.zBridge, false);
		if( closestBelow2Idx!=-1 )
			interp_SmoothPointwise_BetweenTwoConts(closestBelow2Idx, closestBelowIdx);
	}
}

//############################################################


//----------------------------------------------------------------------------
//
//          CONTOUR EVENT FUNCTIONS:
//
//----------------------------------------------------------------------------


//------------------------
//-- Adds a new INTERPOLATED contour to the specified object

bool addInterpolatedContToObj( const InterpolationData &interpolation_data,
															 Iobj *obj, Icont *cont, int interpolated) {
	//## CHECK CONTOUR IS WITHIN TOMOGRRAM'S Z RANGE:

	int zVal = getZ(cont);
	if( zVal < 0 || zVal >= interpolation_data.zsize ) {
		//wprint( "Contour '%d' is outside Z range {%d-%d}\n", zVal, minZ, maxZ );
		return false;
	}

	//## MAKE INTERPOLATED AND ADD:

	Icont *newCont = imodContourDup( cont );    // malloc new contour and don't delele it
	setInterpolated( newCont, interpolated );   // make contour interpolated

	int numConts = csize(obj);
	//undoContourAdditionCO( interpolation_data.view, numConts );             // REGISTER UNDO
	int newContPos = imodObjectAddContour( obj, newCont );
	free(newCont);

	if( newContPos == -1 || newContPos != numConts ) {
		wprint( "ERROR - addInterpolatedContToObj\n" );
		return false;
	}							// FIGURE THIS OUT LATER

	return true;
}

//------------------------
//-- Removes all contours in the object which have their delete flag set to 1

int removeAllDeleteFlaggedContoursFromObj( Iobj *obj )
{
	Icont *cont;
	int numRemoved = 0;
	for( int c=csize(obj)-1; c>=0; c-- )
	{
		cont = getCont(obj, c);
		if( isDeleteFlag( cont ) && isInterpolated( cont ) )
		{
			//undoContourRemovalCO( interpolation_data.view, c );              // REGISTER UNDO
			imodObjectRemoveContour( obj, c );
			numRemoved++;
		}
	}
	return numRemoved;
}





//############################################################



//----------------------------------------------------------------------------
//
//          SIMPLE FUNCTIONS:
//
//----------------------------------------------------------------------------

/*

//---------------------------------
//-- Returns a pointer to the currently selected object.

Iobj *getCurrObj()
{
	Imod *imod = ivwGetModel(interpolation_data.view);
	return ( imodObjectGet(imod) );
}


//---------------------------------
//-- Returns a pointer to the currently selected contour.

Icont *getCurrCont()
{
	Imod *imod = ivwGetModel(interpolation_data.view);
	return ( imodContourGet(imod) );
}

//---------------------------------
//-- Returns a pointer to the currently selected point.

Ipoint *getCurrPt()
{
	Imod *imod = ivwGetModel(interpolation_data.view);
	return ( imodPointGet(imod) );
}


//---------------------------------
//-- Returns true if the object is valid

bool isCurrObjValid()
{
	Iobj *obj = getCurrObj();
	return (obj!=NULL);
}


//---------------------------------
//-- Returns true is a valid contour is selected.

bool isCurrContValid()
{
	return ( isContValid( getCurrCont() ) );
}

//---------------------------------
//-- Returns true is a valid point is selected.

bool isCurrPtValid()
{
	Ipoint *pt = getCurrPt();
	return (pt!=NULL);
}

*/


//----------------------------------------------------------------------------
//
//					POINT INTERPOLATION METHODS:
//
//----------------------------------------------------------------------------


//------------------------
//-- Takes four points for the catmull-rom equation, and through an
//-- iterative high/low guess process tries to resolve the FRACTION
//-- of the distance between p1 and p2 which has the desired Z value.
//--
//-- NOTE: All four points should have their z value in either
//--       ascending or descending order.
//--
//-- NOTE/WARNING:  There is almost certainly a way to compute
//--                this matematically, but I tried using
//--                mathimatica to solve the catmull-rom equation for
//--                "t" and it gave me gibberish

const float  ALLOWABLE_DEV = 0.01;	// allowable deviation of z value from desired value
const int    MAX_GUESSES = 50;			// limits the maximum number of guesses to 50
																		// (although unlikely to ever perform this many)

float calcCardinalSplineFractAtZ( int zVal, Ipoint p0, Ipoint p1,
																						 Ipoint p2,  Ipoint p3 )
{
	if      ( zVal == p1.z )	return 0;		//|- if one of these opposing key points
	else if ( zVal == p2.z )	return 1;		//|  is on the same z value: return it

	if( !( p1.z <= (float)zVal && (float)zVal <= p2.z) ) {		// should never happen
		wprint("ERROR: getPtAndFractWithZValueForCardinalSpline");
		return -1;
	}

	float lowerFract = 0;		// records previous lowest  guess (which value must be > than)
	float upperFract = 1;		// records previous highest guess (which value must be < than)

	float lowerZVal = MIN( p1.z, p2.z );		// Z value of previous lowest  guess
	float upperZVal = MAX( p1.z, p2.z );		// Z value of previous highest guess

	float nextGuessFract;		// the fraction between p1 and p2 (along the spline)
													// of our next guess


	//## USE/RECORD HIGH-LOW GUESSES TO RESOLVE A BEST-GUESS (DOESN'T NEED TO BE EXACT)

	for ( int i=0; i<MAX_GUESSES; i++ )
	{
		nextGuessFract = ((upperFract-lowerFract)/2.0) + lowerFract;
		// sets next guess halfway between the highest and lowest previous safe guesses

		float guessZVal = getValCardinalSpline( nextGuessFract, p0.z, p1.z, p2.z, p3.z,
																				TENSILE_FRACT );
		// determines z value at current guessed fraction

		if (((zVal-ALLOWABLE_DEV)<=guessZVal)       // if estimate is within allowable
				&& (guessZVal<=(zVal+ALLOWABLE_DEV)))   //   deviation of desired z value:
		{
			//cout << " NUM GUESSES:" << i << endl;	//%%%%
			return nextGuessFract;                          // return that fraction
		}
		else if ( guessZVal < zVal )		// else if guess is too low: update lowest guess
		{
			lowerFract = nextGuessFract;
			lowerZVal  = guessZVal;
		}
		else                            // else (if guess is too high):  update highest guess
		{
			upperFract = nextGuessFract;
			upperZVal  = guessZVal;
		}
	}

	return nextGuessFract;
}





//------------------------
//-- Takes four points for the catmull-rom equation, and through an
//-- iterative high/low guess process tries to resolve the FRACTION of all
//-- (integer/slice) z values between p1 and p2 (in ascending order)
//-- and returns this as a vector of floats.
//--
//-- Rather than call calcCardinalSplineFractAtZ (which guesses
//-- and returns a single fract value) this function records previously
//-- guessed values so that the same guess doesn't need to be made twice.
//-- This reduces the number of times getValCardinalSpline is called by about 40%.

std::vector<float> calcCardinalSplineFractsEachSlice( Ipoint p0, Ipoint p1,
																								 Ipoint p2,  Ipoint p3 )
{
	std::vector<float> fracts;

	if ( ABS(p2.z-p1.z) <=1 )   //|- if points p1 and p2 are on the same slice
		return fracts;            //|  or only one apart: return empty set

	int lowerZ = (int)MIN(p1.z,p2.z);				// keeps lowest  Z value of either p1 or p2
	int upperZ = (int)MAX(p1.z,p2.z);				// keeps highest Z value of either p1 or p2

	float lowerZVal = MIN(p1.z,p2.z);		// records Z value at the previous lowest  guess
	float upperZVal = MAX(p1.z,p2.z);		// records Z value at the previous highest guess

	float nextGuessFract;			// the fraction of the way along the path
														// between p1 and p2 of our next guess

	std::vector<IdxToSort> determinedZVals;
	// used to store a list of previous guesses so that the same
	// values won't be fed into getValCardinalSpline twice
	determinedZVals.push_back(IdxToSort(0,lowerZVal,0));
	determinedZVals.push_back(IdxToSort(1,upperZVal,1));
	// we already know these values
	// NOTE: idx = NOT USED, float1 = ZVal of previous guesses,
	//        float2 = actual guess (fract between p1 and p2) )

	int prevLowerBoundInDeterminedList = 0;     // remembers where in the determinedZVals
																							// list the last lowest bound was found

	//## FOR EACH Z VALUE BETWEEN P1 AND P2: USE AND STORE HIGH-LOW GUESSES
	//## TO RESOLVE WHERE THAT Z VALUE LIES

	for ( int targetZ = (lowerZ+1); (targetZ < upperZ); targetZ++ )
		// for all previously recorded (sorted) z values:
	{
		float lowerFract = 0;			// previous lowest  guess (which value must be > than)
		float upperFract = 1;			// previous highest guess (which value must be < than)

		lowerZVal = lowerZ;
		upperZVal = upperZ;

		//## SEACH (SORTED) PREVIOUS GUESS LIST TO FIND CLOSEST GUESS ABOVE
		//## AND BELOW NEW TARGET Z VALUE

		for ( int d=prevLowerBoundInDeterminedList; d<int(determinedZVals.size())-1; d++ )
		{
			if( isBetweenAsc( determinedZVals[d].float1,
												(float)targetZ,
												determinedZVals[d+1].float1 ) )
				// if our desired z value is between these guesses:
				// set starting lowest and highest safe values
			{
				lowerZVal = determinedZVals[d].float1;				//|- (float1 is used to store
				upperZVal = determinedZVals[d+1].float1;			//|  the ZVal of previous guesses)

				lowerFract = determinedZVals[d].float2;				//|- (float2 is used to store
				upperFract = determinedZVals[d+1].float2;			//|  the guesses themselves)

				prevLowerBoundInDeterminedList = d;		// record where in the determinedZVals
																							// value was found, since next z value
																							// will be after this
				break;
			}
		}


		//## USE HIGH-LOW GUESSES TO RESOLVE A BEST-GUESS (DOESN'T NEED TO BE EXACT)

		for ( int i=0; i<MAX_GUESSES; i++ )
		{
			nextGuessFract = ((upperFract-lowerFract)/2.0) + lowerFract;
			// sets next guess as halfway between highest and lowest previous safe guesses

			float guessZVal = getValCardinalSpline( nextGuessFract, p0.z, p1.z, p2.z, p3.z,
																					TENSILE_FRACT );
			// determines z value at current guessed fraction

			if ( isBetweenAsc((targetZ-ALLOWABLE_DEV),guessZVal,(targetZ+ALLOWABLE_DEV)) )
			{		// if estimate is within allowable deviation of desired z value:
					//cout << " NUM GUESSES:" << i << endl;	//%%%%
				break;        // return that fraction
			}
			else if ( guessZVal < targetZ )		// else if guess is too low:
			{
				lowerFract = nextGuessFract;      // update lowest guess
				lowerZVal  = guessZVal;
			}
			else                              // else (if guess is too high):
			{
				upperFract = nextGuessFract;      // update highest guess
				upperZVal  = guessZVal;
				// and add to recorded guesses:
				determinedZVals.push_back( IdxToSort((int)guessZVal,guessZVal,nextGuessFract) );
				// NOTE:  (float1 is used to store the ZVal of previous guesses)
				//        (float2 is used to store the guesses themselves)
			}
		}
		fracts.push_back( nextGuessFract );		// add this best guess to the list of fractions
		determinedZVals = vector_sort( determinedZVals );
	}

	return fracts;
}


