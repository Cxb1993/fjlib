#ifndef FJLIB_CHART_H
#define FJLIB_CHART_H

#include "fjlib_range.h"
#include "qplot_series.h"
#include <boost/shared_ptr.hpp>
#include <list>

namespace qplot {
// really want property feature inside c++
// manage lists of series 
// has zoom build in
// use series to add remove series shared_ptr
// 3.15
// add screen support
struct TQPChart {
	typedef TQPSeries 					series_t;
	typedef boost::shared_ptr<series_t>	series_ptr;
	typedef std::list<series_ptr>		pseries_list;
	
	// properties, i am lazy for the set get crap
	pseries_list
				series;				// series list
	size_t		width,height;		// screen size
//	fjlib::range2_f
//				range;				// series x,y range
	fjlib::float_t
				xy_ratio;			// dx/dy ratio
									// <0 disable ratio
	fjlib::float_t
				margin;				// spacing

	/// default constructor
	TQPChart(): xy_ratio(-1), margin(0),
				width(1), height(1) {}
	/// give you the shared ptr to a new created series
	series_t*	AddNewSeries() {
		series_ptr p(new series_t);
		series.push_back(p);
		return p.get();
	}
	void		SetScreenSize(size_t w, size_t h) {
		width=w; height=h;
	}

	void		TryAdjustRange(const fjlib::range2_f& _range,
					fjlib::float_t& dx, fjlib::float_t& dy) const;
	/// zoom to new range, changes dx or dy range
	fjlib::range2_f
				ZoomXY(const fjlib::range2_f& _range,
					 	fjlib::xyz_direction xy);
	/// zoom to new range, expand range or shink range
	fjlib::range2_f
				ZoomSE(const fjlib::range2_f& _range,
						bool shrink);
	/// zoom in or out
	fjlib::range2_f
				ZoomIn(const fjlib::range2_f& _range,
						fjlib::float_t fraction);
	/// adjust the range based on xy_ratio,
	/// if auto_xy is chosen, xyz_direction returns the final
	/// direction it changes, if auto_xy not chosen,
	///	adjust according to xy direction.
/*
*/
	/// disable xy_ratio
	void		DisableRatio() { xy_ratio=-1; }
};	// end of chart

void TQPChart::TryAdjustRange(const fjlib::range2_f& _range,
					fjlib::float_t& dx, fjlib::float_t& dy) const
{
	if (xy_ratio<=0) { dx=0; dy=0; return; }
	fjlib::float_t ux=width/_range.x.length(),
					uy=height/_range.y.length();
	dy=height/(ux/xy_ratio);
	dx=width/(uy*xy_ratio);
}

fjlib::range2_f TQPChart::ZoomXY(const fjlib::range2_f& _range,
								fjlib::xyz_direction xy)
{
	if (xy_ratio<=0) return _range;
	fjlib::range2_f rg=_range;
	fjlib::float_t dx,dy;
	TryAdjustRange(_range,dx,dy);
	if (xy==fjlib::cdX) 
		rg.x.expand(dx);
	else 
		rg.y.expand(dy);
	return rg;
}

fjlib::range2_f TQPChart::ZoomSE(const fjlib::range2_f& _range,
								bool shrink)
{
	if (xy_ratio<=0) return _range;
	fjlib::range2_f rg=_range;
	fjlib::float_t dx,dy;
	TryAdjustRange(_range,dx,dy);
	fjlib::float_t ddx=dx-_range.x.length();
	bool cx=(shrink == (ddx<=0));
	if (cx)
		rg.x.expand(dx);
	else 
		rg.y.expand(dy);
	return rg;
}

fjlib::range2_f TQPChart::ZoomIn(const fjlib::range2_f& _range,
								fjlib::float_t fraction)
{
	fjlib::range2_f rg=_range;
	fjlib::float_t dx=_range.x.length()*fraction;
	fjlib::float_t dy=_range.y.length()*fraction;
	rg.x.expand(dx);
	rg.y.expand(dy);
	return rg;
}

}	// end of namespace
#endif
