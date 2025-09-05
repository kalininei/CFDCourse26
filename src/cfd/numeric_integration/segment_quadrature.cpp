#include "segment_quadrature.hpp"

using namespace cfd;

const Quadrature* cfd::quadrature_segment_gauss1(){
	static Quadrature quad(
		{
			Point{0}
		},
		{
			2.0
		});
	return &quad;
}

const Quadrature* cfd::quadrature_segment_gauss2(){
	static Quadrature quad(
		{
			Point{-0.5773502691896257},
			Point{0.5773502691896257}
		},
		{
			1.0,
			1.0
		});

	return &quad;
}

const Quadrature* cfd::quadrature_segment_gauss3(){
	static Quadrature quad(
		{
			Point{-0.7745966692414834},
			Point{0},
			Point{0.7745966692414834}
		},
		{
			0.5555555555555556,
			0.8888888888888888,
			0.5555555555555556
		});

	return &quad;
}

const Quadrature* cfd::quadrature_segment_gauss4(){
	static Quadrature quad(
		{
			Point{-0.8611363115940526},
			Point{-0.3399810435848563},
			Point{0.3399810435848563},
			Point{0.8611363115940526}
		},
		{
			0.3478548451374538,
			0.6521451548625461,
			0.6521451548625461,
			0.3478548451374538
		});
	return &quad;
}

const Quadrature* cfd::quadrature_segment_gauss5(){
	static Quadrature quad(
		{
			Point{0.0000000000000000},
			Point{-0.5384693101056831},
			Point{0.5384693101056831},
			Point{-0.9061798459386640},
			Point{0.9061798459386640}
		},
		{
			0.5688888888888889,
			0.4786286704993665,
			0.4786286704993665,
			0.2369268850561891,
			0.2369268850561891
		});
	return &quad;
}


const Quadrature* cfd::quadrature_segment_gauss6(){
	static Quadrature quad(
		{
			Point{0.6612093864662645},
			Point{-0.6612093864662645},
			Point{-0.2386191860831969},
			Point{0.2386191860831969},
			Point{-0.9324695142031521},
			Point{0.9324695142031521}
		},
		{
			0.3607615730481386,
			0.3607615730481386,
			0.4679139345726910,
			0.4679139345726910,
			0.1713244923791704,
			0.1713244923791704
		});
	return &quad;
}
