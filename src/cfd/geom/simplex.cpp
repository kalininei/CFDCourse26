#include "simplex.hpp"

using namespace cfd;

double cfd::triangle_area(Point p0, Point p1, Point p2){
	double x1 = p1.x() - p0.x();
	double y1 = p1.y() - p0.y();
	double x2 = p2.x() - p0.x();
	double y2 = p2.y() - p0.y();

	return 0.5*(x1*y2 - x2*y1);
}
