//#include "fjlib_string.h"
#include "qplot_chart.h"
#include "fjlib_io.h"

int main()
{
	using namespace fjlib;
	range_f x(2,6),y(1,3);
	cout << x.length() << endl;
	range2_f rg(x,y);
	cout << rg.x.length() << endl;

	qplot::TQPChart c;
	c.SetScreenSize(400,320);
	c.xy_ratio=1;
	
	cout << rg << endl << endl;
	cout << rg.within(4,2) << endl;
	cout << c.ZoomXY(rg,fjlib::cdX) << endl;
	cout << c.ZoomXY(rg,fjlib::cdY) << endl;
	cout << c.ZoomSE(rg,true) << endl;
	cout << c.ZoomSE(rg,false) << endl;
	cout << c.ZoomIn(rg,0.5) << endl;
//	cout << c.AdjustRange(range2_f(4,5,4,5)) << endl;
	
//	qplot::TQPChart::series_ptr s=c.AddNewSeries();
//	quickload("in.txt",s->data);
//	cout << c.series.size() << endl;
//	cout << s->data.mat << endl;
/*
	point2_f p(2,3);
	cout << p << endl;

	range2_f rg(x,y);
	cout << rg << endl;
*/
}
