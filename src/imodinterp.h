#pragma once
#include <cstddef>
#include <vector>
#include <iobj.h>
#include <icont.h>
#include <ipoint.h>
#include "an_cont.h"
enum class intmodes : std::ptrdiff_t			    { INT_NO_INTERPOLATE, INT_LINEAR, INT_SPHERICAL, INT_SMOOTH,
																								INT_SMOOTH_PTS, INT_TUBULAR, INT_DEW_DROP, INT_NUM_INT_MODES };

enum class tilingmethod		  { TM_AUTO, TM_CONSERVE_LENGTH, TM_MIN_SA, TM_FEATURE_RECOG };
enum class surfacemethod		{ SR_AUTO, SR_TOUCHING, SR_CENTER_OVERLAP, SR_WITHIN_MIN_DIST,
															SR_FRACTION_OVERLAP, SR_MBR_TOUCH, SR_USER_DEFINED };
enum class branchingmethod	{ BR_BRANCHING_OFF, BR_MERGE_CONTS, BR_BRIDGE_GAPS, BR_PT_BELOW };

enum class ptmethod         { PT_FOUR_PTS, PT_CONVEX_PTS };

constexpr float TENSILE_FRACT = 1.0f;

constexpr int NUM_SAVED_VALS = 13;

struct ContInfo {
	int  idx;             // index of the contour within the object
	int  z;               // the z value of the contour
	bool loaded;          // wether the minimum bounding box has been calculated

	Ipoint ll;            // lower left  point of the contour's minimum bounding box (mbr)
	Ipoint ur;            // upper right point of the contour's minimum bounding box
	Ipoint centerPt;      // point in the center of the MBR

	int checked;          // can be used to see which contour is checked

	void resetAll();
};

struct ZList {
	std::vector<int> idxs;     // a vector of contour indexes
};

struct InterpolationData {
	//ImodView    *view;			// NOTE: members are in /include/imodP.h
	//Interpolator *window;

	intmodes interType;          // the type of interpolation the user is using (see enum "intmodes")
	tilingmethod tilingMethod;       // the type of tiling method being used (see enum "tilingmethod")
	branchingmethod branchingMethod;		// the type of branching method being used to generted interpolated contours (see enum "branchingmethod")
	surfacemethod surfResolveMethod;	// the type of surface resolution method being used to determine which contours are part of the same surface (see enum "surfacemethod")
	ptmethod ptResolveMethod;    // the type of point connecting method used when tilingMethod is set to TM_CONSERVE_LENGTH (see enum "ptmethod")

	int zBridge;                    // number of slices over which to connect contours for interpolation.
	float interSepDistBetweenConts;	// the X-Y distance contours are allowed to be apart, but still considered part of the same surface (when surfResolveMethod == SR_WITHIN_MIN_DIST).
	float interFractOverlap;        // percentage of smaller contours which must overlap the bigger one to be considered in the same same surface (when surfResolveMethod == SR_FRACTION_OVERLAP).
	bool interLineTracker;          // if 1, then line tracker is applied to generated/interpolated contours ( %%%%%%%%%%% NOT YET IMPLEMENTED %%%%%%%%%%% )
	int minHoleSize;                // minimum hole size to find when [h] is pressed
	int maxGapSize;                 // maximum distance used in "edit_findNextIsolatedContour"
	bool hideSurfSettings;					// if true: will hide surface settings area

	//** OTHER:

	bool deselectAfterEnter;    // wether the current contour is deselected after [Enter]
	int selectedAction;         // the last selected action under "More Actions"

	int contIdxLastInterp;			// the contour index of the last contour interpolated
	int objIdxLastInterp;				// the object  index of the last contour interpolated

	bool initialized;           // is set to true after values have been set
	int xsize, ysize, zsize;    // size of the image / tomogram
};

class imod_interpolation {
public:
	imod_interpolation();
	void resetAll();
	void regenerateContInfoVector( int _minZLimit=0, int _maxZLimit=MAX_INT, int startIdx=0 );
	void regenerateZTables( );
	void regenerateConnections( );

	//** USER CALLED METHODS:
	bool performInterpolation(Iobj * const out_obj, Imod * const imod);
	void performInterpolationOnCont( Iobj* _obj, int _contIdx, const intmodes interpolationType, const float modScaleZ = 1.f);
	int  deleteImmediatelyAdjacentInterpolatedConts( Iobj* _obj, int _contId );
	int  findMiddleNextLargeInterpolatedSpan( Iobj *_obj, int _contIdx, int minHoleSize );
	//** SURFACE RESOLVING METHODS:

	bool contoursSameSurf( int c1Idx, int c2Idx );
	bool contoursSameSurfSlow( Icont *cont1, Icont *cont2 );

	std::vector<int> findIdxsNearestKeyConts( int bcIdx, bool above, int maxZDist, int minZDist=1  );
	int findIdxClosestKeyContInList( int bcIdx, std::vector<int> conts );

	int findIdxNearestKeyCont( int cbIdx, int maxDist, bool above );
	std::vector<int> findIdxsAllKeyContsNoBranching( int idxBaseCont, int maxDist );
	int deleteAllSameSurfInterpContsInZRange( int minZ, int maxZ, int cbIdx );
	int deleteAllSameSurfInterpContsEitherSideOfCont( int cbIdx );
	int deleteInterpolatedContsBetweenKeyContsAboveAndBelow();

	//** POINT FITTING METHODS:

	void findCoorrespondingPts( Icont *cont1, Icont *cont2, int *idxCont1, int *idxCont2 );
	void findCoorrespondingPts_FourPtMBR( Icont *cont1, Icont *cont2, int *idxCont1, int *idxCont2 );
	void findCoorrespondingPts_AllConvexPts( Icont *cont1, Icont *cont2, int *idxCont1, int *idxCont2 );      // NEW

	void breakBaseContIntoTwoForBranching( Icont* baseCont, Icont* contBranch1, Icont* contBranch2, Icont *newBranch1, Icont *newBranch2 );

	void modifyKeyContsForInterp( int cbIdx, int ctIdx, Icont *contLNew, Icont *contUNew, bool closeContours, bool findCoorrespondingPtsAndMakeClockwise );
	void modifyKeyContsForInterp_LineConservation( Icont *contKeyL, Icont *contKeyU, Icont *contL, Icont *contU, bool closeContours, bool findCoorrespondingPtsAndMakeClockwise );
	void modifyKeyContsForInterp_MinimumAreaCost ( int cbIdx, int ctIdx, Icont *contLNew, Icont *contUNew );


	//** INTERPOLATION METHODS:

	void interp_SmoothCrude_BetweenConts( int c0, int c1, int c2, int c3 );
	void interp_SmoothPointwise_BetweenTwoConts( int clIdx, int cuIdx );
	std::vector<IcontPtr> getLinearInterpConts( int c1Idx, int c2Idx );
	int interp_Linear_BetweenTwoConts( int clIdx, int cuIdx );

	int  mergeAllTouchingConts( std::vector<IcontPtr> conts );
	void addAllInterpolatedConstAndMerge( std::vector< std::vector<IcontPtr> > newConts );
	void interp_Linear_BetweenContsMerge( int cbIdx, std::vector<int> branchContIdx );
	int interp_Linear( int baseContIdx, int maxDist );
	void interp_Spherical_OnSingleCont( int cMiddleIdx, const float modScaleZ );
	void interp_Spherical( int baseContIdx, int maxDist, const float modScaleZ );
	void interp_SmoothCrude( int baseContIdx, int maxDist );
	void interp_SmoothCrudeAll( int baseContIdx, int maxDist );
	void interp_SmoothPointwise( int baseContIdx, int maxDist );


	//** SIMPLE INLINE ACCESSORS:

	inline int numKeyContsAtZ(int z)         {  return static_cast<int>(ztableKey[std::size_t(z)].idxs.size()); }
	inline int numIntContsAtZ(int z)         {  return static_cast<int>(ztableInt[std::size_t(z)].idxs.size()); }

	inline int idxKeyContZ(int z, int i)     {  return ztableKey[std::size_t(z)].idxs[std::size_t(i)]; }
	inline int idxIntContZ(int z, int i)     {  return ztableInt[std::size_t(z)].idxs[std::size_t(i)]; }

	inline bool isKeyContAbove()             { return ( closestAboveIdx != -1 );	}
	inline bool isKeyContBelow()             { return ( closestBelowIdx != -1 );	}

	inline Icont*  getC(const int idx) const								 { return ( getCont(obj, idx) ); }
	inline int     getCZ(const std::size_t idx) const				 { return ( conti[idx].z ); }
	inline Ipoint* getUR(const std::size_t idx)							 { return ( &conti[idx].ur ); }
	inline Ipoint* getLL(const std::size_t idx)						   { return ( &conti[idx].ll ); }
	inline Ipoint* getCenterPt(const std::size_t idx)		     { return ( &conti[idx].centerPt ); }
	inline bool    getClosed(int idx) const									 { return ( isContClosed(obj, getC(idx)) ); }
public:
	//## DATA:

		Iobj  * obj;           // pointer to the OBJECT  we are dealing with
		int sContIdx;         // index of the contour we must interpolate (within *obj)
		int sContZ;           // the slice with *cont is on (the starting z value)

		std::vector<ContInfo> conti;     // is populated with contour info (such as MBR) for each
																// contour matching it's order in *obj

		int minZLimit;
		int maxZLimit;

		std::vector<ZList> ztableKey;    // list of indexes of key contours on each slice
		std::vector<ZList> ztableInt;    // list of indexes of interpolated contours on each slice

		int closestAboveIdx;  // the NEAREST same-surface contour ABOVE our given contour
		int closestBelowIdx;	// the NEAREST same-surface contour BELOW our given contour

		std::vector<int> aboveIdxs;		// list of NEAREST contours ABOVE given contour with
															//   same z val - (if > 1 then given contour branches)
		std::vector<int> belowIdxs;		// list of NEAREST contours BELOW given contour with
															//   same z val - (if > 1 then given contour branches)
		InterpolationData interpolation_data;
};

//-------------------------------
//## SMALL FUNCTIONS:

Iobj *getCurrObj();
Icont *getCurrCont();
Ipoint *getCurrPt();
bool isCurrObjValid();
bool isCurrContValid();
bool isCurrPtValid();

//-------------------------------
//## CONTOUR EVENT FUNCTIONS:

bool addInterpolatedContToObj(const InterpolationData &interpolation_data, Iobj *obj, Icont *cont, int interpolated=1);
int removeAllDeleteFlaggedContoursFromObj( Iobj *obj );


//-------------------------------
//## POINT INTERPOLATION FUNCTIONS:

float calcCardinalSplineFractAtZ( int zVal, Ipoint p0, Ipoint p1, Ipoint p2,  Ipoint p3 );
std::vector<float> calcCardinalSplineFractsEachSlice( Ipoint p0, Ipoint p1, Ipoint p2,  Ipoint p3 );

