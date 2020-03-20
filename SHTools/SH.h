
#pragma once

#include "SHEval.h"
#include <cmath>

#define PI 3.1415927f

static constexpr int SH_ORDER = 3;
static constexpr int SH_COUNT = SH_ORDER * SH_ORDER;

enum FILTER_TYPE
{
	FILTER_DISABLE,
	FILTER_GAUSSIAN,
	FILTER_HANNING,
	FILTER_LANCZOS
};

class SHCoefficients
{
public:
	float& operator[](int i)
	{
		return val[i];
	}

	const float& operator[](int i) const
	{
		return val[i];
	}

	float* getPtr()
	{
		return val;
	}

private:
	float val[SH_COUNT] = { 0 };
};

class SHColor
{
public:
	SHCoefficients& operator[](int i)
	{
		return sh[i];
	}

	const SHCoefficients& operator[](int i) const
	{
		return sh[i];
	}

private:
	SHCoefficients sh[3];
};

void projectOnSH(float x, float y, float z, SHCoefficients& sh)
{
	//static constexpr float Y00 = 0.28209479177387814f;	// sqrt(1/4pi)
	//static constexpr float Y10 = 0.4886025119029199f;	// sqrt(3/4pi)
	//static constexpr float Y21 = 1.0925484305920792f;	// sqrt(15/4pi)
	//static constexpr float Y20 = 0.31539156525252005f;	// sqrt(5/16pi)
	//static constexpr float Y22 = 0.5462742152960396f;	// sqrt(15/16pi)

	//sh[0] = Y00;								// L = 0, M = 0
	//sh[1] = -Y10 * y;							// L = 1, M = -1
	//sh[2] = Y10 * z;							// L = 1, M = 0
	//sh[3] = -Y10 * x;							// L = 1, M = 1
	//sh[4] = Y21 * x * y;						// L = 2, M = -2
	//sh[5] = -Y21 * y * z;						// L = 2, M = -1
	//sh[6] = Y20 * (3.0f * z * z - 1.0f);		// L = 2, M = 0
	//sh[7] = -Y21 * x * z;						// L = 2, M = 1
	//sh[8] = Y22 * (x * x - y * y);				// L = 2, M = 2

	float* pSH = sh.getPtr();
	switch (SH_ORDER)
	{
	case 3:
		SHEval3(x, y, z, pSH);
		break;
	case 4:
		SHEval4(x, y, z, pSH);
		break;
	case 5:
		SHEval5(x, y, z, pSH);
		break;
	case 6:
		SHEval6(x, y, z, pSH);
		break;
	case 7:
		SHEval7(x, y, z, pSH);
		break;
	case 8:
		SHEval8(x, y, z, pSH);
		break;
	case 9:
		SHEval9(x, y, z, pSH);
		break;
	case 10:
		SHEval10(x, y, z, pSH);
		break;

	default:
		break;
	}
}

void reconstructFromSH(float x, float y, float z, const SHColor& shColor, float color[3])
{
	SHCoefficients shDirection;
	projectOnSH(x, y, z, shDirection);

	color[0] = color[1] = color[2] = 0.0f;
	for (int c = 0; c < 3; ++c)
	{
		for (int i = 0; i < SH_COUNT; ++i)
		{
			color[c] += shColor[c][i] * shDirection[i];
		}
	}
}

void filter(int l, float a, SHCoefficients& sh)
{
	int offset = l * l;
	for (int m = -l; m <= l; ++m, ++offset)
	{
		sh[offset] *= a;
	}
}

void gauss(int w, SHCoefficients& sh)
{
	for (int l = 0; l < SH_ORDER; ++l)
	{
		filter(l, std::exp(-std::pow(PI * float(l) / float(w), 2.0f) / 2.0f), sh);
	}
}

void hanning(int w, SHCoefficients& sh)
{
	for (int l = 0; l < SH_ORDER; ++l)
	{
		if (l > w)
			filter(l, 0, sh);
		else
			filter(l, (std::cos(PI * float(l) / float(w)) + 1.0f) * 0.5f, sh);
	}
}

void lanczos(int w, SHCoefficients& sh)
{
	for (int l = 0; l < SH_ORDER; ++l)
	{
		if (l == 0)
			filter(l, 1, sh);
		else
			filter(l, std::sin(PI * float(l) / float(w)) / (PI * float(l) / float(w)), sh);
	}
}
